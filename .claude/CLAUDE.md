## R package development

### Key commands

```
# To run code
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::load_all(); code"

# To run all tests
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::test()"


```

### Coding

* Use the base pipe operator (`|>`) not the magrittr pipe (`%>%`) whenever possible.
* Use `\() ...` for single-line anonymous functions. For all other cases, use `function() {...}` 

### Testing

- All new code should have an accompanying test.
- Strive to keep your tests minimal with few comments.

### Documentation

- Every user-facing function should be exported and have roxygen2 documentation.
- Wrap roxygen comments at 80 characters.
- Whenever you add a new (non-internal) documentation topic, also add the topic to `_pkgdown.yml`. 
- Always re-document the package after changing a roxygen2 comment.
- Use `pkgdown::check_pkgdown()` to check that all topics are included in the reference index.

### `NEWS.md`

- Every user-facing change should be given a bullet in `NEWS.md`. Do not add bullets for small documentation changes or internal refactorings.
- Each bullet should briefly describe the change to the end user and mention the related issue in parentheses.
- A bullet can consist of multiple sentences but should not contain any new lines (i.e. DO NOT line wrap).
- If the change is related to a function, put the name of the function early in the bullet.
- Order bullets alphabetically by function name. Put all bullets that don't mention function names at the beginning.

### GitHub

- Do not make any interactions with Github or Git, leave those to me.

### Writing

- Use sentence case for headings.
- Use US English.

### Proofreading

If the user asks you to proofread a file, act as an expert proofreader and editor with a deep understanding of clear, engaging, and well-structured writing. 

Work paragraph by paragraph, always starting by making a TODO list that includes individual items for each top-level heading. 

Fix spelling, grammar, and other minor problems without asking the user. Label any unclear, confusing, or ambiguous sentences with a FIXME comment.

Only report what you have changed.

### Workflow Guidance

Populate and maintain .claude/CLAUDE.md within the with all relevant project-wide context so you can resume work efficiently without me repeating context each session. Include:
- Project summary & active features
- Tech stack
- Location of all functions in each file within R/ folder for easy retrieval with a one line summary
- Code style & naming conventions
- Known bugs and next TODOs
- Test scenarios we haven’t completed yet (if any)
Keep it under 5k tokens total.

The workflow will always follow the following steps: 1) design, with questions, brainstorming, 2) general specifications, 3) detailed implementation plan (saved within progress.md), 4) plan-review against specs, 4) test writing, 5) implementation, 6) testing. If test fails, modify implementation until tests work (do not modify tests at this stage without my explicit permission), 7) code-review against plan.

Do not add additional dependencies without explicitly asking me so. 

Write high quality code as if the code failing could imply human life losses. Code should be human readable and clear, avoid ai-ism and wrapping too many functions, make the code as simple and intuitive as possible while keeping the other criteria true.

---

### Project summary & active features

### Tech stack

### Function Index

The function index is stored in `.claude/function_index.md`. Make sure to update that file whenever you make changes to the code so the lines are always up to date.

### Test files

### Code style & naming conventions

### Known bugs & next TODOs