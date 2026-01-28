# Submodule forks for muhammadhasyim

**fairchem** already uses your repo: `https://github.com/muhammadhasyim/fairchem.git` — no fork needed.

**cace** and **LES-BEC** point to BingqingCheng’s repos. To use your own forks under `muhammadhasyim`:

## 1. Create forks on GitHub

1. **cace**  
   - Open: https://github.com/BingqingCheng/cace  
   - Click **Fork** (top right) → choose owner **muhammadhasyim** → Create fork.

2. **LES-BEC**  
   - Open: https://github.com/BingqingCheng/LES-BEC  
   - Click **Fork** → choose owner **muhammadhasyim** → Create fork.

## 2. Push local submodule commits to your forks

From the **openmm** repo root:

```bash
# cace: add your fork, push main
cd cace
git remote add myfork https://github.com/muhammadhasyim/cace.git 2>/dev/null || true
git push myfork main

# LES-BEC: add your fork, push main
cd ../LES-BEC
git remote add myfork https://github.com/muhammadhasyim/LES-BEC.git 2>/dev/null || true
git push myfork main

cd ..
```

## 3. Switch submodule URLs to your forks (optional)

So future `git submodule update --init` and clones use your forks:

```bash
git submodule set-url cace https://github.com/muhammadhasyim/cace.git
git submodule set-url LES-BEC https://github.com/muhammadhasyim/LES-BEC.git
git submodule sync
```

Then commit the updated `.gitmodules`:

```bash
git add .gitmodules && git commit -m "Point cace and LES-BEC submodules to muhammadhasyim forks"
```
