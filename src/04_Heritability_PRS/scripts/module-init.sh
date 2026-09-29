# Sourced only inside the command launchers; no changes to the RStudio process.
if ! type module >/dev/null 2>&1; then
  if [[ -r /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
  else
    echo 'Biowulf module initialization is unavailable. Use tool_mode: path for a prepared installation.' >&2
    exit 2
  fi
fi
