# Conda Forge releases for OpenMM

## Final releases

Create a new version tag on `openmm/openmm` and publish a new GitHub release. The Conda Forge bots would auto submit a PR with the new version within a couple hours _if_ we were using the `source.url` field pointing to a GH release. Since we build off the full repo, I am not sure the automation will work, so we need to switch to manual mode.

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. On your `conda-forge/openmm-feedstock` **fork**, create a new branch off most recent upstream branch, be it `master` or `rc`.
3. Edit `meta.yaml`:
   - [ ] Update `package.version`.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.git_rev` points to the release tag you want to publish. This could be a git commit too, but a tag is preferred.
4. Commit and push to **your fork**. Do NOT push to `upstream`.
5. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targeting `conda-forge/openmm-feedstock`'s `master`, from **your fork**.
6. Review the checklist and open the PR.
7. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
8. Wait for all green, reviews and then merge. Always make sure you are merging to `master`.
9. Once merged, check the CI status on the `master` branch.

   - It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
   - If a CI provider was not triggered (for whatever reason), `master` might need a little _push_ (no pun intended). An empty commit to `master` will do:

   ```
   # make sure conda-forge/openmm-feedstock is configured as `upstream`
   git remote -v
   git checkout master
   git fetch upstream master
   git merge upstream/master
   git commit --allow-empty -m "Trigger CI"
   git push upstream master
   ```

## Release candidates

> Technically, once you publish an RC tag on GitHub Releases, the bots will pick it up, but we haven't tested this yet.

Manual instructions:

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. Create a new branch from the most recent upstream branch, be it `master` or `rc`.
3. Edit `meta.yaml`:
   - [ ] Update `package.version`. It should be the new version number plus `rcX`, `X` being a number. Check [CFEP05](https://github.com/conda-forge/cfep/blob/master/cfep-05.md) in case of doubt.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.git_rev` points to the release tag you want to publish. This could be a git commit too, but a tag is preferred.
4. Edit `conda_build_config.yaml`. This a KEY difference: RC releases are published to a different label!

   - [ ] `channel_targets` should be set to `- conda-forge openmm_rc`:

   ```yaml
   channel_targets:
     - conda-forge openmm_rc
   ```

5. Commit and push to **your fork**. Do NOT push to `upstream`.
6. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targeting `conda-forge/openmm-feedstock`'s `rc`, from **your fork**. Again, we are **targeting** the `rc` branch, NOT master. This is a KEY difference. RC candidates stay on `rc`.
7. Review the checklist and open the PR.
8. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
9. Wait for all green, reviews and then merge. Always make sure you are merging to `rc`.
10. Once merged, check the CI status on the `rc` branch.

- It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
- If a CI provider was not triggered (for whatever reason), `rc` might need a little _push_ (no pun intended). An empty commit to `rc` will do:

```
# make sure conda-forge/openmm-feedstock is configured as `upstream`
git remote -v
git checkout rc
git fetch upstream rc
git merge upstream/rc
git commit --allow-empty -m "Trigger CI"
git push upstream rc
```

## Development Releases

These releases are done "on demand" either if users request them or if we feel there are some important or useful code changes we want people to use or try out.
An advantage these releases have over the nightly omnia builds is that they are ABI compatible with the conda-forge ecosystem.
They are uploaded to conda-forge in a way that prevents users from accidentally installing them.
See [CFEP-05](https://github.com/conda-forge/cfep/blob/main/cfep-05.md) and [PR# 74](https://github.com/conda-forge/openmm-feedstock/pull/74) for more details.

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
1. Create a new branch from the upstream `dev` branch.
1. Edit `meta.yaml`"
   - [ ] Update `package.version`. If there has not been a new openmm release since the last dev release, bump the dev build number. Example: If current release of openmm is `7.7.0` and the package version is currently `7.7.0dev0` change it to `7.7.0dev1`. Another example, if the current release of openmm is `7.7.1` and the package version is currently `7.7.0dev2` then change it to `7.7.1dev0`.
   - [ ] Update the `source.git_rev` to the current HEAD commit of openmm.
1. Commit and push to **your fork**. Do NOT push to `upstream`.
1. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targeting `conda-forge/openmm-feedstock`'s `dev`, from **your fork**. Again, we are **targeting** the `dev` branch, NOT master. This is a KEY difference.
1. Review the checklist and open the PR.
1. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
1. Wait for all green, reviews and then merge. Always make sure you are merging to `dev`.
1. Once merged, check the CI status on the `dev` branch.

- It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
- If a CI provider was not triggered (for whatever reason), `dev` might need a little _push_ (no pun intended). An empty commit to `dev` will do:

```
# make sure conda-forge/openmm-feedstock is configured as `upstream`
git remote -v
git checkout rc
git fetch upstream dev
git merge upstream/dev
git commit --allow-empty -m "Trigger CI"
git push upstream dev
```

## Hosted Environments

Anaconda.org supports hosting environments (as defined in a yaml) in the Anaconda cloud.
This makes it easier for end users to create a specified conda environment as they do not need a local copy of the yaml file.
This is especially useful when distributing an environment that uses lots of channels and/or labels.

To upload en environment

1. Create the yaml file
1. Test that the yaml file works e.g. `mamba env create --file openmm-8-beta.yaml`
1. Upload the environment e.g. `anaconda upload -u openmm hosted-envs/openmm-8-beta.yaml`

Now users can create an environment with `mamba env create openmm/openmm-8-beta`
