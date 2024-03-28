# Contributing Guide
 Contributions are welcome and appreciated. You can contribute by reporting a bug, fixing a bug, implementing a new feature, or writing documentation.

## Report Bugs
Report a bug by opening a new issue. In your bug report, you should include:
* A quick summary and/or background
* Operating system name and version
* Steps to reproduce
    * Be specific
    * Give sample code if possible
* Expected behavior
* Actual behavior
* Notes
    * Why you think this might be happening
    * Things you tried that did not work

## Fix Bugs
Fix any bug in the GitLab project tagged with "bug" or "help wanted."

## Implement Features
Any issue in the GitLab project tagged with "enhancement" and that is unassigned is open to anyone to implement.

## Write Documentation
Update the documents, classes, or libraries with up-to-date or more detailed information.

## Getting Started
1. Clone the GitLab project.
2. Create a branch for local development that corresponds to a specific feature or a significant bug fix. For minor fixes
, no need to make a new branch.
3. Install requirements with `make init` (only need to do this once)
5. Contribute your changes
6. Compile Cython files with `make compile`
    * You need to do this anytime a Cython file (*.pyx) is created or modified.
7. Run tests with `make test`
8. Ensure your code is formatted correctly by running `make reformat`

The final steps vary depending on your intent.
#### Commit to Python3 Branch
1. Add to the [changelog](../CHANGELOG.MD) and follow semantic versioning.
2. Commit your changes. See [commit guidelines](#commits) for more.
3. Submit a merge request and use the "MERGE_TEMPLATE" template.

#### Commit to Local Branch
If the tests pass and your code is formatted correctly, commit your changes. See the guidelines below for more.
  
### Commits
Make sure that your commit message is
 specific. You can find commit advice [here](https://thoughtbot.com/blog/5-useful-tips-for-a-better-commit-message
 ) and [here](https://dev.to/jacobherrington/how-to-write-useful-commit-messages-my-commit-message-template-20n9
 ). Additionally, you can use the commit template below by adding it to your .gitconfig:
 ```
########50 characters############################
Subject
########72 characters##################################################
Task
# Problem, Task, Reason for Commit
Change(s):
# Solution or List of Changes
Note
# Special instructions, testing steps, rake, etc

```
