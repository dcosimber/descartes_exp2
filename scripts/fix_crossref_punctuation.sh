#!/usr/bin/env bash
set -euo pipefail

if [[ ! -d docs ]]; then
  exit 0
fi

find docs -type f -name '*.html' -print0 | xargs -0 perl -0pi -e '
  s/(Tabla&nbsp;(?:<span>)?[0-9.]+(?:<\/span>)?)[:.]{1,2}\s+/$1. /g;
  s/(Figura&nbsp;(?:<span>)?[0-9.]+(?:<\/span>)?)[:.]{1,2}\s+/$1. /g;
'
