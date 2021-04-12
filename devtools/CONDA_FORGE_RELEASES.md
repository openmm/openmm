# Conda Forge releases for OpenMM

## Final releases

Create a new version tag on `openmm/openmm` and publish a new GitHub release. The Conda Forge bots will auto submit a PR with the new version within a couple hours. Go to step 8 below.

> If the automated PR doesn't happen for any reason (wait 24h just in case), follow the manual instructions below.

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. On your `conda-forge/openmm-feedstock` **fork**, create a new branch off latest upstream's `master`.
3. Edit `meta.yaml`:
   - [ ] Update `package.version`.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.url` is correct.
   - [ ] Update `source.sha256` with the output of `openssl sha256 downloaded_file_from_source_url.tar.gz`
4. Commit and push to **your fork**. Do NOT push to `upstream`.
5. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targetting `conda-forge/openmm-feedstock`'s `master`, from **your fork**.
6. Review the checklist and open the PR.
7. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
8. Wait for all green, reviews and then merge. Always make sure you are merging to `master`.
9. Once merged, check the CI status on the `master` branch. It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).

## Release candidates

> Technically, once you publish an RC tag on GitHub Releases, the bots will pick it up, but we haven't tested this yet.

Manual instructions:

1. If you haven't yet, fork `conda-forge/openmm-feedstock`.
2. If needed, on your `conda-forge/openmm-feedstock` **fork**, sync the `rc` branch so it contains the latest changes in `master`.
3. Create a new branch off latest upstream's `rc`. This is a KEY difference with final releases. RC candidates stay on `rc`.
4. Edit `meta.yaml`:
   - [ ] Update `package.version`. It should be the new version number plus `rcX`, `X` being a number. Check [CFEP05](https://github.com/conda-forge/cfep/blob/master/cfep-05.md) in case of doubt.
   - [ ] Reset `build.number` to 0.
   - [ ] Make sure `source.url` is correct.
   - [ ] Update `source.sha256` with the output of `openssl sha256 downloaded_file_from_source_url.tar.gz`
5. Edit `conda_build_config.yaml`:

   - [ ] `channel_targets` should be set to `- conda-forge openmm_rc`:

   ```yaml
   channel_targets:
     - conda-forge openmm_rc
   ```

6. Commit and push to **your fork**. Do NOT push to `upstream`.
7. Open a new PR on `conda-forge/openmm-feedstock`. Make sure you are targetting `conda-forge/openmm-feedstock`'s `rc`, from **your fork**. Again, we are **targetting** the `rc` branch, NOT master.
8. Review the checklist and open the PR.
9. In the opened PR, post a comment with `@conda-forge-admin, please rerender`.
10. Wait for all green, reviews and then merge. Always make sure you are merging to `rc`.
11. Once merged, check the CI status on the `rc` branch. It should be green. If it's red, a network error might have happened and you need to _re-run_ the failing job (a link will appear next to it, if you click on the Details menu).
