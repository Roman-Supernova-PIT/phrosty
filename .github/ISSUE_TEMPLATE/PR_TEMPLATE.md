---
name: "Pull Request for phrosty"
about: Create a pull request to submit new or updated code to the repository
title: "[PR]: "
labels: possible solution
assignees: ''
---
## What type of PR is this? (check all applicable)

- [ ] Refactor
- [ ] Feature
- [ ] Bug Fix
- [ ] Optimization
- [ ] Documentation Update
- [ ] Have you followed the guidelines in our Contributing document?
- [ ] Have you checked to ensure there aren't other open [Pull Requests](../../../pulls) for the same update/change?

<!-- describe the changes included with this PR here -->
This PR addresses ...

<!-- If this PR closes a GitHub issue, reference it here by its number -->
Closes #


<!-- if you can't perform these tasks due to permissions, please ask a maintainer to do them -->
## Tasks
- [ ] **request a review from someone specific**, to avoid making the maintainers review every PR
- [ ] Does this PR change user-facing code / API? (if not, label with `no-changelog-entry-needed`)
  - [ ] write news fragment(s) in `changes/`: `echo "changed something" > changes/<PR#>.<changetype>.rst` (see expansion below for change types) 
  - [ ] update or add relevant tests
  - [ ] update relevant docstrings and / or `docs/` page
    - [ ] Do truth files need to be updated?


<details><summary>news fragment change types...</summary>

- ``changes/<PR#>.general.rst``: infrastructure or miscellaneous change
- ``changes/<PR#>.docs.rst``: updates for documentation
</details>