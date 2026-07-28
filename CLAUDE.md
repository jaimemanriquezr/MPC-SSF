# CLAUDE.md

## Plans and decisions

Three record-keeping conventions. **This flow is mandatory for non-trivial work.**

### `.claude/plans/` — one file per refactor or feature

Named `YYYY-MM-DD-slug.md`. Contains:

- **Goal** — what "done" means.
- **Steps** — each with the verification check that proves it worked.
- **Files touched.**

Written **before** implementing, not after. If the approach changes mid-flight,
update the plan rather than letting it drift out of date.

### `.claude/decisions/` — one file per significant design decision

Named `YYYY-MM-DD-slug.md`. Contains:

- **The decision.**
- **Why** — the domain or engineering reasoning.
- **Alternatives rejected, and why.**

This is the searchable "why". Cite evidence with paths (`julia/analysis/...`,
`src/simulate.jl:385`). When a decision rests on a measurement, the numbers go in
a durable file that the decision links to — **never only in the conversation**.
A decision that cites a number with no file behind it is not finished.

### `.claude/JOURNAL.md` — chronological development log

One dated section per session: what was done, what was learned, what's next.

Complementary to `decisions/`: the **journal records what happened**, the
**decision files record why**. A session that produced a decision gets both — a
journal line noting it, and the decision file carrying the reasoning.

### The rule

For any non-trivial work:

1. **Plan before implementing** — write the plan file first.
2. **Decision file whenever an alternative was genuinely considered** — if you
   weighed two approaches and picked one, that reasoning gets written down.
3. **Journal entry before ending a session.**

Trivial work (typo fixes, one-line corrections, mechanical edits) is exempt. When
unsure whether something qualifies, write the record — the cost of an unnecessary
file is far lower than the cost of an unrecoverable "why".
