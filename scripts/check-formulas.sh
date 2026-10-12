#!/usr/bin/env bash
# Cross-check the formula tags.
#
#   Every Rust function that implements a formula carries a doc line  [formula: <kebab-id>].
#   Each such id must be either
#     * on a Lean theorem (or lemma) in lean/proofs  -> status "proved", or
#     * listed in lean/proofs/CITED.md               -> status "cited".
#   Anything else is "MISSING" and makes the script fail.
#
#   The script also fails if `sorry`, `admit` or an `axiom` declaration appears in the Lean
#   sources of lean/proofs.
#
# Usage:  scripts/check-formulas.sh [--build]
#   --build   additionally run `lake build` in lean/proofs and check that the axiom audit
#             (DdgProofs/Audit.lean) reports only propext, Classical.choice and Quot.sound.
#
# Works with no Rust sources present (the table then lists the Lean and cited ids only).

set -u
LC_ALL=C
export LC_ALL

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
proofs="$root/lean/proofs"
cited="$proofs/CITED.md"
crates="$root/crates"
do_build=0
for arg in "$@"; do
  case "$arg" in
    --build) do_build=1 ;;
    -h|--help) sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
fail=0

# ---------------------------------------------------------------------------------------------
# 1. Tags on Lean declarations:  id <TAB> kind <TAB> name <TAB> file:line
#    A tag belongs to the first declaration that follows it.
# ---------------------------------------------------------------------------------------------
: > "$tmp/lean.tsv"
if [ -d "$proofs" ]; then
  find "$proofs" -name '*.lean' -not -path '*/.lake/*' -print | sort | while IFS= read -r f; do
    rel="${f#"$root"/}"
    awk -v file="$rel" '
      function flush(kind, name, line,   i) {
        for (i = 1; i <= n; i++) printf "%s\t%s\t%s\t%s:%d\n", ids[i], kind, name, file, line
        n = 0
      }
      {
        s = $0
        while (match(s, /\[formula: [a-z0-9-]+\]/)) {
          id = substr(s, RSTART + 10, RLENGTH - 11)
          ids[++n] = id
          s = substr(s, RSTART + RLENGTH)
        }
        t = $0
        sub(/^[ \t]+/, "", t)
        while (match(t, /^(private|protected|noncomputable|nonrec|unsafe|partial|@\[[^]]*\])[ \t]+/))
          t = substr(t, RLENGTH + 1)
        if (match(t, /^(theorem|lemma|def|abbrev|structure|inductive|instance|class|example)[ \t(:]/)) {
          kind = t; sub(/[ \t(:].*$/, "", kind)
          name = substr(t, length(kind) + 1); sub(/^[ \t]+/, "", name); sub(/[ \t({\[:].*$/, "", name)
          if (kind == "example" || name == "") name = "(anonymous)"
          if (n > 0) flush(kind, name, NR)
        }
      }
      END { if (n > 0) flush("none", "(no declaration follows)", NR) }
    ' "$f" >> "$tmp/lean.tsv"
  done
fi

# ---------------------------------------------------------------------------------------------
# 2. Cited ids:  id <TAB> CITED.md:line
# ---------------------------------------------------------------------------------------------
: > "$tmp/cited.tsv"
if [ -f "$cited" ]; then
  grep -n -o '\[formula: [a-z0-9-]*\]' "$cited" 2>/dev/null \
    | sed -E 's/^([0-9]+):\[formula: ([a-z0-9-]+)\]$/\2\tlean\/proofs\/CITED.md:\1/' \
    | sort -u -t "$(printf '\t')" -k1,1 > "$tmp/cited.tsv"
fi

# ---------------------------------------------------------------------------------------------
# 3. Rust tags:  id <TAB> file:line
# ---------------------------------------------------------------------------------------------
: > "$tmp/rust.tsv"
if [ -d "$crates" ]; then
  find "$crates" -name '*.rs' -not -path '*/target*/*' -print0 2>/dev/null \
    | xargs -0 -r grep -n -o '\[formula: [a-z0-9-]*\]' /dev/null 2>/dev/null \
    | sed -E "s|^$root/||" \
    | sed -E 's/^(.*):([0-9]+):\[formula: ([a-z0-9-]+)\]$/\3\t\1:\2/' > "$tmp/rust.tsv"
fi

# ---------------------------------------------------------------------------------------------
# 4. The table
# ---------------------------------------------------------------------------------------------
{
  cut -f1 "$tmp/rust.tsv"
  cut -f1 "$tmp/lean.tsv"
  cut -f1 "$tmp/cited.tsv"
} | sort -u > "$tmp/ids.txt"

n_proved=0; n_cited=0; n_missing=0; n_rust=0; n_unused=0
printf '%-38s %-9s %-5s %s\n' "id" "status" "rust" "location"
printf '%-38s %-9s %-5s %s\n' "--------------------------------------" "---------" "-----" "--------"
while IFS= read -r id; do
  [ -n "$id" ] || continue
  in_rust=no
  if grep -q "^$id	" "$tmp/rust.tsv"; then in_rust=yes; n_rust=$((n_rust + 1)); fi
  thm="$(awk -F '\t' -v id="$id" '$1 == id && ($2 == "theorem" || $2 == "lemma") { print $4 "  " $3 }' "$tmp/lean.tsv")"
  other="$(awk -F '\t' -v id="$id" '$1 == id && !($2 == "theorem" || $2 == "lemma") { print $4 "  " $3 " (" $2 ", not a theorem)" }' "$tmp/lean.tsv")"
  cit="$(awk -F '\t' -v id="$id" '$1 == id { print $2 }' "$tmp/cited.tsv")"
  if [ -n "$thm" ]; then
    status=proved; n_proved=$((n_proved + 1)); loc="$thm"
  elif [ -n "$cit" ]; then
    status=cited; n_cited=$((n_cited + 1)); loc="$cit"
  else
    status=MISSING; loc="${other:-(no Lean theorem, not in CITED.md)}"
    if [ "$in_rust" = yes ]; then n_missing=$((n_missing + 1)); fail=1; else status=untagged; fi
  fi
  [ "$in_rust" = no ] && n_unused=$((n_unused + 1))
  first=1
  printf '%s\n' "$loc" | while IFS= read -r line; do
    if [ "$first" = 1 ]; then
      printf '%-38s %-9s %-5s %s\n' "$id" "$status" "$in_rust" "$line"; first=0
    else
      printf '%-38s %-9s %-5s %s\n' "" "" "" "$line"
    fi
  done
  if [ "$in_rust" = yes ]; then
    awk -F '\t' -v id="$id" '$1 == id { printf "%-38s %-9s %-5s   used by %s\n", "", "", "", $2 }' "$tmp/rust.tsv"
  fi
done < "$tmp/ids.txt"

echo
echo "ids: $(wc -l < "$tmp/ids.txt" | tr -d ' ')   proved: $n_proved   cited: $n_cited   used by Rust: $n_rust   missing: $n_missing   not (yet) used by Rust: $n_unused"

if [ "$n_missing" -gt 0 ]; then
  echo "FAIL: $n_missing Rust formula id(s) have neither a Lean theorem nor an entry in CITED.md" >&2
fi

# An id that is both proved and cited is suspicious (stale CITED.md entry).
while IFS= read -r id; do
  [ -n "$id" ] || continue
  if awk -F '\t' -v id="$id" '$1 == id && ($2 == "theorem" || $2 == "lemma") { found = 1 } END { exit !found }' "$tmp/lean.tsv" \
     && grep -q "^$id	" "$tmp/cited.tsv"; then
    echo "note: '$id' is on a Lean theorem and also listed in CITED.md" >&2
  fi
done < "$tmp/ids.txt"

# Tags that are not attached to any declaration.
if grep -q '	none	' "$tmp/lean.tsv"; then
  echo "FAIL: formula tag(s) in lean/proofs not followed by a declaration:" >&2
  grep '	none	' "$tmp/lean.tsv" | awk -F '\t' '{ print "  " $1 "  " $4 }' >&2
  fail=1
fi

# ---------------------------------------------------------------------------------------------
# 5. No sorry, no admit, no axiom declarations in lean/proofs
# ---------------------------------------------------------------------------------------------
bad=0
if [ -d "$proofs" ]; then
  find "$proofs" -name '*.lean' -not -path '*/.lake/*' -print0 \
    | xargs -0 -r grep -n -E '(^|[^[:alnum:]_.])(sorry|admit)([^[:alnum:]_]|$)' /dev/null > "$tmp/sorry.txt" 2>/dev/null
  find "$proofs" -name '*.lean' -not -path '*/.lake/*' -print0 \
    | xargs -0 -r grep -n -E '^[[:space:]]*(@\[[^]]*\][[:space:]]*)?((private|protected|noncomputable|unsafe)[[:space:]]+)*axiom[[:space:]]' /dev/null > "$tmp/axiom.txt" 2>/dev/null
  if [ -s "$tmp/sorry.txt" ]; then
    echo "FAIL: 'sorry'/'admit' found in lean/proofs:" >&2; sed "s|^$root/|  |" "$tmp/sorry.txt" >&2; bad=1
  fi
  if [ -s "$tmp/axiom.txt" ]; then
    echo "FAIL: axiom declaration found in lean/proofs:" >&2; sed "s|^$root/|  |" "$tmp/axiom.txt" >&2; bad=1
  fi
fi
if [ "$bad" = 0 ]; then
  echo "no sorry, no admit, no axiom declarations in lean/proofs"
else
  fail=1
fi

# ---------------------------------------------------------------------------------------------
# 6. Optional: build and audit the axioms actually used
# ---------------------------------------------------------------------------------------------
if [ "$do_build" = 1 ]; then
  mkdir -p "$root/runs_tmp"
  log="$root/runs_tmp/check-formulas-build.log"
  echo "building lean/proofs (log: runs_tmp/check-formulas-build.log) ..."
  if (cd "$proofs" && lake build DdgProofs > "$log" 2>&1 && lake env lean DdgProofs/Audit.lean >> "$log" 2>&1); then
    if grep -q "sorryAx" "$log"; then
      echo "FAIL: the axiom audit mentions sorryAx" >&2; fail=1
    fi
    extra="$(tr '\n' ' ' < "$log" | grep -o 'depends on axioms: \[[^]]*\]' | sed -E 's/depends on axioms: \[(.*)\]/\1/' | tr ',' '\n' | sed 's/^ *//; s/ *$//' | sort -u | grep -v -x -E 'propext|Classical\.choice|Quot\.sound|' || true)"
    if [ -n "$extra" ]; then
      echo "FAIL: non-standard axioms in use: $extra" >&2; fail=1
    else
      echo "build ok; $(grep -c 'depends on axioms\|does not depend on any axioms' "$log") audited theorems use only propext, Classical.choice, Quot.sound"
    fi
  else
    echo "FAIL: lake build failed, see $log" >&2; fail=1
  fi
fi

if [ "$fail" = 0 ]; then echo "OK"; else echo "FAILED" >&2; fi
exit "$fail"
