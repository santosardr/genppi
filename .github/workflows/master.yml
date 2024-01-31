name: Lisp CI

on:
  push:
    branches:
      - master
  pull_request:
    branches:
      - master

jobs:
  build:
    runs-on: macos-latest

    steps:
    - name: Checkout repository
      uses: actions/checkout@v2

    - name: Install SBCL
      run: brew install sbcl

    - name: Setup Quicklisp
      run: |
        curl -O https://beta.quicklisp.org/quicklisp.lisp
        sbcl --non-interactive --load quicklisp.lisp --eval '(quicklisp-quickstart:install)'

    - name: Cache Quicklisp packages
      uses: actions/cache@v2
      with:
        path: ~/quicklisp/
        key: ${{ runner.os }}-quicklisp-${{ hashFiles('**/*.asd', '**/*.lisp') }}
        restore-keys: |
          ${{ runner.os }}-quicklisp-

    - name: Load dependencies and local projects
      run: |
        sbcl --non-interactive \
             --load ~/quicklisp/setup.lisp \
             --eval '(push (merge-pathnames "src/" (user-homedir-pathname)) ql:*local-project-directories*)' \
             --eval '(ql:quickload :lparallel)' \
             --eval '(ql:quickload :cl-random-forest)' \
             --eval '(ql:quickload :cl-store)' \
             --eval '(ql:quickload :features)'

    - name: Build and save the executable
      run: |
        sbcl --noinform --dynamic-space-size 8192 --control-stack-size 20 \
             --load src/genppi.lisp \
             --eval "(save-lisp-and-die \"binaries/genppi8g-Mac.x\" :executable t :save-runtime-options t :toplevel 'main)"

    - name: Upload the macOS executable as an artifact
      uses: actions/upload-artifact@v2
      with:
        name: genppi8g-Mac
        path: binaries/genppi8g-Mac.x