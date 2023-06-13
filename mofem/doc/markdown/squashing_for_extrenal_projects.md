# How to: Squashing for an external PR {#SquashFAQ}
## 1. What is squashing?
Squashing combines multiple commits into one, in order to keep a large project's commit history compact.
## 2. When to squash?
a. Before making a PR  
b. After making additional changes as per request of developers
## 3. How to:
1. Checkout the branch which you (will) PR from:   
    ```
    git checkout feature_1
    ```
2. Merge all changes into the feature branch (local):  
    ```
    git merge feature_1_dev
    ```
3. Rebase with respect to the upstream master of the third party repo:   
    ```
    git rebase -i remotes/upstream/master
    ```
4. A list of all commits will appear in reverse chronological order (\# comments out a line):   
    ```
    pick bcf4f3fa58 commit oldest
    pick 7bfec3e3dc commit 1
    pick d65272ef0d commit 2
    pick c4d48c1856 commit 3
    pick d1bd7824bc commit 4
    pick e43f01d7ea commit latest
    ```
5. One commit, usually the original (oldest) commit, has to be left as a pick and all succeeding/preceding commits must be squashed into it, using the ```squash``` or ```s``` key word as follows:
    ```
    pick bcf4f3fa58 commit oldest
    s 7bfec3e3dc commit 1
    s d65272ef0d commit 2
    s c4d48c1856 commit 3
    s d1bd7824bc commit 4
    s e43f01d7ea commit latest
    ```
6. All commit messages will appear together, you must then edit them accordingly. The resulting text will be the squashed commit message.  
    * If any merge errors occur the rebase will be paused and you will be able to resolve them. Once you have selected the correct changes you need to add the problem file using:  
        ```
        git add merge_error.cpp
        ```
        Then the rebase process must be continued with:
        ```
        git rebase --continue
        ```
7. The new commit needs to be pushed to the remote branch using:
    ```
    git push -f
    ```
    The ```-f``` means force and is needed to replace the remote commit history with the new local commit history. More info on this can be found at: [https://stackoverflow.com/questions/39399804/updates-were-rejected-because-the-tip-of-your-current-branch-is-behind-its-remot]
    