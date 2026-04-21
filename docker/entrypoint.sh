#!/bin/sh
set -eu
cd "${BAYES_CRICKET_ROOT:-/app}"
# Docker often runs without a TTY, so R block-buffers stdout and logs look "empty".
# script(1) runs the command with a pseudo-tty; -e propagates the child's exit status;
# -f flushes output; typescript is discarded to /dev/null.
exec script -qefc "Rscript R/run_example.R" /dev/null
