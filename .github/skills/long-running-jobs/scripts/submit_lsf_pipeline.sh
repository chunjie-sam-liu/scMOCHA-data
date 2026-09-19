#!/usr/bin/env bash
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-16
# @DESCRIPTION: Submit a multi-stage LSF pipeline as one dependency graph so the
#   whole chain queues in a single call and runs unattended overnight. The DAG is
#   read from the Mermaid block in the campaign's PROGRESS.md, so the graph the
#   user reads is the graph that was submitted. Fan-out and fan-in work today.
#   Also refreshes a submitted graph (--status) and reprints the Mermaid block.
#   Usage: submit_lsf_pipeline.sh [options] --graph <PROGRESS.md>
#          submit_lsf_pipeline.sh [options] <stage.lsf>...
#          submit_lsf_pipeline.sh --status <state.tsv>
# @VERSION: v0.2.0

set -euo pipefail

PROG=${0##*/}

die() {
  printf '%s: error: %s\n' "$PROG" "$*" >&2
  exit 1
}

warn() {
  printf '%s: warning: %s\n' "$PROG" "$*" >&2
}

if ((BASH_VERSINFO[0] < 4)); then
  die "bash 4+ required (associative arrays); found ${BASH_VERSION}"
fi

usage() {
  cat <<'EOF'
Usage:
  submit_lsf_pipeline.sh [options] --graph <PROGRESS.md>
  submit_lsf_pipeline.sh [options] <stage.lsf> [<stage.lsf> ...]
  submit_lsf_pipeline.sh [options] --spec <pipeline.tsv>
  submit_lsf_pipeline.sh --status <state.tsv> [--no-mermaid]

Submits every stage at once and lets LSF hold each one PEND until its
prerequisites finish successfully.

Sources for the DAG, in order of preference:
  --graph   the first ```mermaid fence in a markdown file, normally the
            campaign PROGRESS.md. The graph the user reads is the graph that
            gets submitted, so there is nothing to keep in sync.
  positional .lsf files, which form a linear chain
  --spec    a standalone 3-column file, for a pipeline that has no PROGRESS.md
            yet or for testing

Graph contract (flowchart TD or LR):
  node id      the stage id. [A-Za-z0-9_] only
  node label   fields separated by <br/>. The field holding a *.lsf or
               *.sbatch path makes it a stage; the rest is display text
  edge         parent --> child, with or without an edge label
  ext_<jobid>  a job already in the queue
  any other node (no script path) is a display-only artifact box and is
  dropped together with its edges

  stage_id["export<br/>track/run_export.lsf [1-2]<br/>done 2/2"]:::done

Options:
  --graph FILE       Read the DAG from the first mermaid fence in FILE.
  --spec FILE        DAG spec. Whitespace-separated columns:
                       id  script.lsf  deps  [extra bsub args]
                     deps is a comma-separated list of ids, or - for none.
  --after IDS        Comma-separated LSF job IDs already in the queue. Every
                     root stage waits on all of them.
  --after-array      Treat every --after ID as a job array (skip the bjobs
                     probe). Use when the probe cannot reach mbatchd.
  --from IDS         Resume: submit only these stages and their descendants.
  --skip IDS         Treat these stages as already finished; drop them and
                     rewire their children onto their dependencies.
  --state FILE       Where to write the submitted graph.
                     Default: logs/lsf-pipeline/<name>-<timestamp>.tsv
  --status FILE      Read-only. Re-read a state file, query LSF, print the
                     current table, a refreshed Mermaid block, and the exact
                     recovery commands. Submits nothing.
  -n, --dry-run      Print the submission plan and exit. Nothing is submitted.
      --no-mermaid   Suppress the Mermaid block.
  -h, --help         This message.

Dependency expressions (see the lsf-resource-planner skill):
  array parent    numdone(<jobid>,*)   every element DONE, i.e. all exited 0
  single parent   done(<jobid>)        that job DONE

  NEVER done(<jobid>[*]). In LSF that is a one-to-one element-wise dependency
  requiring both arrays to be the same size, not "all elements succeeded".

This script writes nothing back into the markdown. It prints the updated block;
pasting it into PROGRESS.md is a normal edit. It never runs bkill either: when
a stage fails, its downstream jobs sit PEND on a condition that can never be
satisfied, and --status prints the bkill and --from commands for you to review.
EOF
}

# ---------------------------------------------------------------- graph model

# A dependency token is either an internal stage id, or @<jobid> for a job that
# is already in the queue.
declare -a ORDER=()
declare -A SCRIPT=() DEPS=() EXTRA=() JOBID=() JNAME=() ISARRAY=() CHILDREN=()
declare -A EXT_ISARRAY=() EXT_CLASS=() CLASS=() NOTE=()

add_stage() {
  local id=$1 script=$2 deps=$3 extra=${4:-}
  [[ -z ${SCRIPT[$id]:-} ]] || die "duplicate stage id '${id}'"
  [[ $id =~ ^[A-Za-z0-9_]+$ ]] ||
    die "stage id '${id}' must be [A-Za-z0-9_] only, so it is a valid Mermaid node id"
  ORDER+=("$id")
  SCRIPT[$id]=$script
  DEPS[$id]=$deps
  EXTRA[$id]=$extra
}

# The trailing newline is load-bearing: `while read` drops a final line that has
# none, which silently turns a one-item dependency list into no dependency.
split_csv() {
  printf '%s\n' "${1:-}" | tr ',' '\n' | sed '/^[[:space:]]*$/d;/^-$/d'
}

kind_word() {
  if (($1)); then printf 'array'; else printf 'single job'; fi
}

# ------------------------------------------------------------- argument parse

spec_file=""
graph_file=""
after_raw=""
after_array=0
from_raw=""
skip_raw=""
state_file=""
status_file=""
dry_run=0
want_mermaid=1
declare -a lsf_files=()

while (($# > 0)); do
  case $1 in
    --spec)
      spec_file=${2:?--spec needs a file}
      shift 2
      ;;
    --graph)
      graph_file=${2:?--graph needs a markdown file}
      shift 2
      ;;
    --after)
      after_raw=${2:?--after needs one or more job IDs}
      shift 2
      ;;
    --after-array)
      after_array=1
      shift
      ;;
    --from)
      from_raw=${2:?--from needs one or more stage ids}
      shift 2
      ;;
    --skip)
      skip_raw=${2:?--skip needs one or more stage ids}
      shift 2
      ;;
    --state)
      state_file=${2:?--state needs a file}
      shift 2
      ;;
    --status)
      status_file=${2:?--status needs a state file}
      shift 2
      ;;
    -n | --dry-run)
      dry_run=1
      shift
      ;;
    --no-mermaid)
      want_mermaid=0
      shift
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      die "unknown option '$1' (--help for usage)"
      ;;
    *)
      lsf_files+=("$1")
      shift
      ;;
  esac
done
lsf_files+=("$@")

# ------------------------------------------------------------------ rendering

# The emitted block is the input format: the node id is the stage id and the
# label carries the script path, so a printed graph can be submitted again.
# The six classDef colors are fixed documentation chrome for PROGRESS.md, not a
# figure palette: copy them verbatim and do not move them into a track color.R.
mermaid_block() {
  local id dep node
  local -A seen_ext=()
  printf '```mermaid\n'
  printf 'flowchart TD\n'
  for id in "${ORDER[@]}"; do
    printf '  %s["%s<br/>%s %s<br/>%s"]:::%s\n' \
      "$id" "$id" "${SCRIPT[$id]}" "${JNAME[$id]:-}" \
      "${NOTE[$id]:-not submitted}" "${CLASS[$id]:-todo}"
  done
  for id in "${ORDER[@]}"; do
    while read -r dep; do
      [[ $dep == @* ]] || continue
      [[ -z ${seen_ext[$dep]:-} ]] || continue
      seen_ext[$dep]=1
      printf '  ext_%s["already queued<br/>%s"]:::%s\n' \
        "${dep#@}" "${dep#@}" "${EXT_CLASS[${dep#@}]:-done}"
    done < <(split_csv "${DEPS[$id]}")
  done
  for id in "${ORDER[@]}"; do
    while read -r dep; do
      if [[ $dep == @* ]]; then
        node="ext_${dep#@}"
      else
        node=$dep
      fi
      printf '  %s --> %s\n' "$node" "$id"
    done < <(split_csv "${DEPS[$id]}")
  done
  cat <<'EOF'

  classDef done fill:#D0E9E6,stroke:#2A9D8F,stroke-width:2px,color:#1B665D
  classDef run  fill:#CCE3F2,stroke:#1982C4,stroke-width:2px,color:#10557F
  classDef pend fill:#FCE6D0,stroke:#F28E2B,stroke-width:2px,color:#9D5C1C
  classDef todo fill:#E3DFDE,stroke:#BAB0AC,stroke-width:1px,color:#79706E
  classDef fail fill:#FADFD9,stroke:#E76F51,stroke-width:3px,color:#B13E06
  classDef dead fill:#DED8E7,stroke:#6A4C93,stroke-width:2px,stroke-dasharray:4 3,color:#453160
EOF
  printf '```\n'
}

# ------------------------------------------------------------------ lsf facts

# LSF purges finished job records after CLEAN_PERIOD, so an empty bjobs answer
# means "record gone", not "not an array". Return 2 so the caller cannot
# conflate the two.
probe_is_array() {
  local jid=$1 name
  name=$(bjobs -noheader -o "job_name" "$jid" 2>/dev/null | head -n1 || true)
  [[ -n $name ]] || return 2
  [[ $name == *\[*\]* ]]
}

jname_of() {
  local id=$1 tok prev="" j=""
  local -a toks=()
  read -r -a toks <<<"${EXTRA[$id]:-}" || true
  for tok in ${toks[@]+"${toks[@]}"}; do
    if [[ $prev == "-J" ]]; then j=$tok; fi
    prev=$tok
  done
  if [[ -z $j ]]; then
    j=$(awk '$1 == "#BSUB" { for (i = 2; i <= NF; i++) if ($i == "-J") { print $(i + 1); exit } }' \
      "${SCRIPT[$id]}")
  fi
  j=${j%\"}
  j=${j#\"}
  j=${j%\'}
  j=${j#\'}
  printf '%s' "$j"
}

dep_term() {
  local jid=$1 is_array=$2
  if ((is_array)); then
    printf 'numdone(%s,*)' "$jid"
  else
    printf 'done(%s)' "$jid"
  fi
}

# bjobs -A columns: JOBID ARRAY_SPEC OWNER NJOBS PEND DONE RUN EXIT ...
array_counts() {
  local jid=$1 line stat
  line=$(bjobs -A -noheader "$jid" 2>/dev/null | head -n1 || true)
  if [[ -n $line ]]; then
    awk '{ print $4, $5, $6, $7, $8 }' <<<"$line"
    return 0
  fi
  stat=$(bjobs -noheader -o "stat" "$jid" 2>/dev/null | head -n1 || true)
  case $stat in
    DONE) printf '1 0 1 0 0\n' ;;
    EXIT) printf '1 0 0 0 1\n' ;;
    RUN) printf '1 0 0 1 0\n' ;;
    PEND | PSUSP | USUSP | SSUSP) printf '1 1 0 0 0\n' ;;
    *) return 1 ;;
  esac
}

# --------------------------------------------------------------- build graph

load_spec() {
  local file=$1 line id script deps rest lineno=0
  while IFS= read -r line || [[ -n $line ]]; do
    lineno=$((lineno + 1))
    line=${line%%#*}
    [[ -n ${line//[[:space:]]/} ]] || continue
    read -r id script deps rest <<<"$line"
    [[ -n $script ]] || die "${file}:${lineno}: need at least 'id script deps'"
    add_stage "$id" "$script" "${deps:--}" "${rest:-}"
  done <"$file"
}

load_linear() {
  local prev="" f id
  for f in "${lsf_files[@]}"; do
    id=${f##*/}
    id=${id%.lsf}
    id=${id#run_}
    add_stage "$id" "$f" "${prev:--}"
    prev=$id
  done
}

# Read the DAG out of the first ```mermaid fence in a markdown file. A node
# whose label carries no *.lsf path is a display-only artifact box and is
# dropped along with its edges, so a graph may carry more than the submission.
load_graph() {
  local file=$1 line bare id label field src dst i
  local in_fence=0
  local -a node_ids=() src_list=() dst_list=()
  local -A gscript=() gdeps=() gkind=()

  while IFS= read -r line || [[ -n $line ]]; do
    bare=${line#"${line%%[![:space:]]*}"}
    if ((in_fence == 0)); then
      if [[ $bare == '```mermaid'* ]]; then in_fence=1; fi
      continue
    fi
    [[ $bare != '```'* ]] || break
    case $bare in
      classDef* | class\ * | %%* | flowchart* | graph\ * | subgraph* | end | style\ * | linkStyle* | click*) continue ;;
    esac
    if [[ $bare =~ ^([A-Za-z0-9_]+)\[\"(.*)\"\] ]]; then
      id=${BASH_REMATCH[1]}
      label=${BASH_REMATCH[2]}
      [[ -z ${gkind[$id]:-} ]] || die "${file}: node '${id}' is defined twice in the graph"
      gkind[$id]=display
      node_ids+=("$id")
      while IFS= read -r field; do
        case $field in
          *.lsf* | *.sbatch*) ;;
          *) continue ;;
        esac
        gscript[$id]=$(awk '{ for (i = 1; i <= NF; i++) if ($i ~ /\.(lsf|sbatch)$/) { print $i; exit } }' <<<"$field")
        [[ -z ${gscript[$id]} ]] || gkind[$id]=stage
        break
      done < <(printf '%s\n' "${label//<br\/>/$'\n'}")
      if [[ ${gkind[$id]} == display && $id =~ ^ext_[0-9]+$ ]]; then gkind[$id]=ext; fi
    elif [[ $bare == *--*\>* ]]; then
      src=$(awk '{ print $1 }' <<<"${bare%%-*}")
      dst=$(awk '{ print $1 }' <<<"${bare##*>}")
      if [[ -n $src && -n $dst ]]; then
        src_list+=("$src")
        dst_list+=("$dst")
      fi
    fi
  done <"$file"

  ((in_fence)) || die "${file}: no \`\`\`mermaid fence found"

  for ((i = 0; i < ${#src_list[@]}; i++)); do
    src=${src_list[$i]}
    dst=${dst_list[$i]}
    [[ ${gkind[$dst]:-} == stage ]] || continue
    case ${gkind[$src]:-} in
      stage) gdeps[$dst]="${gdeps[$dst]:+${gdeps[$dst]},}${src}" ;;
      ext) gdeps[$dst]="${gdeps[$dst]:+${gdeps[$dst]},}@${src#ext_}" ;;
    esac
  done

  for id in "${node_ids[@]}"; do
    [[ ${gkind[$id]} == stage ]] || continue
    add_stage "$id" "${gscript[$id]}" "${gdeps[$id]:--}"
  done
  ((${#ORDER[@]} > 0)) ||
    die "${file}: the first mermaid fence has no node whose label carries a .lsf path"
}

load_state() {
  local file=$1 stage script jname jid deps rest
  while IFS=$'\t' read -r stage script jname jid deps rest; do
    [[ -n $stage && $stage != \#* ]] || continue
    add_stage "$stage" "$script" "${deps:--}"
    JNAME[$stage]=$jname
    if [[ -n $jid && $jid != "-" ]]; then JOBID[$stage]=$jid; fi
  done <"$file"
  ((${#ORDER[@]} > 0)) || die "no stages found in ${file}"
}

# ---------------------------------------------------------------- graph utils

build_children() {
  local id dep
  CHILDREN=()
  for id in "${ORDER[@]}"; do
    while read -r dep; do
      [[ $dep != @* ]] || continue
      [[ -n ${SCRIPT[$dep]:-} ]] || die "stage '${id}' depends on unknown stage '${dep}'"
      CHILDREN[$dep]="${CHILDREN[$dep]:-} $id"
    done < <(split_csv "${DEPS[$id]}")
  done
}

topo_sort() {
  local -A indeg=()
  local -a queue=() out=()
  local id dep cur ch
  for id in "${ORDER[@]}"; do indeg[$id]=0; done
  for id in "${ORDER[@]}"; do
    while read -r dep; do
      [[ $dep != @* ]] || continue
      indeg[$id]=$((indeg[$id] + 1))
    done < <(split_csv "${DEPS[$id]}")
  done
  for id in "${ORDER[@]}"; do
    if ((indeg[$id] == 0)); then queue+=("$id"); fi
  done
  while ((${#queue[@]} > 0)); do
    cur=${queue[0]}
    queue=(${queue[@]+"${queue[@]:1}"})
    out+=("$cur")
    for ch in ${CHILDREN[$cur]:-}; do
      indeg[$ch]=$((indeg[$ch] - 1))
      if ((indeg[$ch] == 0)); then queue+=("$ch"); fi
    done
  done
  ((${#out[@]} == ${#ORDER[@]})) ||
    die "the spec has a cycle: only ${#out[@]} of ${#ORDER[@]} stages are reachable"
  ORDER=("${out[@]}")
}

# Drop a stage and rewire: its children inherit its dependencies.
drop_stage() {
  local gone=$1 id dep merged ch
  local -a kept=()
  for id in "${ORDER[@]}"; do
    [[ $id != "$gone" ]] || continue
    merged=""
    while read -r dep; do
      if [[ $dep == "$gone" ]]; then
        while read -r ch; do
          merged="${merged:+${merged},}${ch}"
        done < <(split_csv "${DEPS[$gone]}")
      else
        merged="${merged:+${merged},}${dep}"
      fi
    done < <(split_csv "${DEPS[$id]}")
    DEPS[$id]=${merged:--}
    kept+=("$id")
  done
  ORDER=(${kept[@]+"${kept[@]}"})
  unset "SCRIPT[$gone]" "DEPS[$gone]" "EXTRA[$gone]"
}

descendants_of() {
  local -A seen=()
  local -a stack=("$@")
  local cur ch
  while ((${#stack[@]} > 0)); do
    cur=${stack[0]}
    stack=(${stack[@]+"${stack[@]:1}"})
    [[ -z ${seen[$cur]:-} ]] || continue
    seen[$cur]=1
    for ch in ${CHILDREN[$cur]:-}; do stack+=("$ch"); done
  done
  if ((${#seen[@]} > 0)); then printf '%s\n' "${!seen[@]}"; fi
}

# ===================================================================== status

if [[ -n $status_file ]]; then
  [[ -r $status_file ]] || die "cannot read state file '${status_file}'"
  load_state "$status_file"
  build_children
  topo_sort

  declare -A NJOBS=() NPEND=() NDONE=() NRUN=() NEXIT=()
  for id in "${ORDER[@]}"; do
    NJOBS[$id]=0 NPEND[$id]=0 NDONE[$id]=0 NRUN[$id]=0 NEXIT[$id]=0
    if [[ -z ${JOBID[$id]:-} ]]; then
      CLASS[$id]=todo
      continue
    fi
    if read -r n p d r e < <(array_counts "${JOBID[$id]}"); then
      NJOBS[$id]=$n NPEND[$id]=$p NDONE[$id]=$d NRUN[$id]=$r NEXIT[$id]=$e
      if ((e > 0)); then
        CLASS[$id]=fail
      elif ((r > 0)); then
        CLASS[$id]=run
      elif ((n > 0 && d == n)); then
        CLASS[$id]=done
      else
        CLASS[$id]=pend
      fi
    else
      CLASS[$id]=todo
      JNAME[$id]="${JNAME[$id]:-}(cleaned)"
    fi
  done

  # A stage still PEND behind a failed or dead ancestor can never be satisfied.
  for id in "${ORDER[@]}"; do
    [[ ${CLASS[$id]} == pend ]] || continue
    while read -r dep; do
      [[ $dep != @* ]] || continue
      case ${CLASS[$dep]:-} in
        fail | dead)
          CLASS[$id]=dead
          break
          ;;
      esac
    done < <(split_csv "${DEPS[$id]}")
  done

  printf '\nPipeline status: %s\n\n' "$status_file"
  printf '%-16s %-12s %-24s %-5s %5s %5s %5s %5s\n' \
    STAGE JOBID NAME STATE DONE RUN PEND EXIT
  for id in "${ORDER[@]}"; do
    printf '%-16s %-12s %-24s %-5s %5s %5s %5s %5s\n' \
      "$id" "${JOBID[$id]:--}" "${JNAME[$id]:--}" "${CLASS[$id]}" \
      "${NDONE[$id]}" "${NRUN[$id]}" "${NPEND[$id]}" "${NEXIT[$id]}"
    NOTE[$id]="${JOBID[$id]:--} ${CLASS[$id]} ${NDONE[$id]}/${NJOBS[$id]}"
  done

  declare -a broken=() orphan=()
  for id in "${ORDER[@]}"; do
    case ${CLASS[$id]} in
      fail) broken+=("$id") ;;
      dead) orphan+=("${JOBID[$id]}") ;;
    esac
  done

  if ((${#broken[@]} > 0)); then
    printf '\nBroken at: %s\n' "${broken[*]}"
    printf 'Read the first failures before changing anything:\n'
    for id in "${broken[@]}"; do
      printf '  grep -l . logs/*/*_%s_*.err 2>/dev/null | head -3 | xargs -r -n1 tail -n 30\n' \
        "${JOBID[$id]}"
    done
    if ((${#orphan[@]} > 0)); then
      printf '\nThese downstream jobs are PEND on a condition that can never be\n'
      printf 'satisfied. Review, then kill them before resubmitting:\n'
      printf '  bkill %s\n' "${orphan[*]}"
    fi
    printf '\nAfter fixing, rerun only the failed indices, then resume the tail:\n'
    printf '  bsub -J "%s[<failed-indices>]" < %s\n' \
      "${JNAME[${broken[0]}]%%[*}" "${SCRIPT[${broken[0]}]}"
    printf '  %s --graph <PROGRESS.md> --after <new-jobid> --from %s\n' \
      "$PROG" "$(printf '%s' "${CHILDREN[${broken[0]}]:-<next-stage>}" | awk '{ print $1 }')"
  fi

  if ((want_mermaid)); then
    printf '\n'
    mermaid_block
  fi
  exit 0
fi

# ===================================================================== submit

if [[ -n $graph_file ]]; then
  [[ -z $spec_file && ${#lsf_files[@]} -eq 0 ]] ||
    die "use exactly one of --graph, --spec, or positional .lsf files"
  [[ -r $graph_file ]] || die "cannot read graph '${graph_file}'"
  load_graph "$graph_file"
elif [[ -n $spec_file ]]; then
  ((${#lsf_files[@]} == 0)) || die "use exactly one of --graph, --spec, or positional .lsf files"
  [[ -r $spec_file ]] || die "cannot read spec '${spec_file}'"
  load_spec "$spec_file"
elif ((${#lsf_files[@]} > 0)); then
  load_linear
else
  usage >&2
  exit 2
fi

build_children
topo_sort

# --skip and --from are both expressed as stage removals, so the DAG stays valid.
declare -a to_drop=()
while read -r id; do
  [[ -n ${SCRIPT[$id]:-} ]] || die "--skip names unknown stage '${id}'"
  to_drop+=("$id")
done < <(split_csv "$skip_raw")

if [[ -n $from_raw ]]; then
  declare -a roots=()
  while read -r id; do
    [[ -n ${SCRIPT[$id]:-} ]] || die "--from names unknown stage '${id}'"
    roots+=("$id")
  done < <(split_csv "$from_raw")
  declare -A keep=()
  while read -r id; do keep[$id]=1; done < <(descendants_of "${roots[@]}")
  for id in "${ORDER[@]}"; do
    [[ -n ${keep[$id]:-} ]] || to_drop+=("$id")
  done
fi

for id in ${to_drop[@]+"${to_drop[@]}"}; do
  [[ -n ${SCRIPT[$id]:-} ]] || continue
  drop_stage "$id"
  build_children
done
((${#ORDER[@]} > 0)) || die "--from / --skip removed every stage"

# Job IDs already in the queue become dependencies of every remaining root.
declare -a after_ids=()
while read -r jid; do
  [[ $jid =~ ^[0-9]+$ ]] || die "--after '${jid}' is not a numeric LSF job ID"
  after_ids+=("$jid")
done < <(split_csv "$after_raw")

# The graph may already name ext_<jobid> nodes; those need the same probe.
declare -a probe_ids=()
declare -A probe_seen=()
for jid in ${after_ids[@]+"${after_ids[@]}"}; do
  if [[ -z ${probe_seen[$jid]:-} ]]; then
    probe_seen[$jid]=1
    probe_ids+=("$jid")
  fi
done
for id in "${ORDER[@]}"; do
  while read -r dep; do
    [[ $dep == @* ]] || continue
    jid=${dep#@}
    if [[ -z ${probe_seen[$jid]:-} ]]; then
      probe_seen[$jid]=1
      probe_ids+=("$jid")
    fi
  done < <(split_csv "${DEPS[$id]}")
done

for jid in ${probe_ids[@]+"${probe_ids[@]}"}; do
  EXT_CLASS[$jid]=done
  if ((after_array)); then
    EXT_ISARRAY[$jid]=1
    continue
  fi
  rc=0
  probe_is_array "$jid" || rc=$?
  case $rc in
    0) EXT_ISARRAY[$jid]=1 ;;
    1) EXT_ISARRAY[$jid]=0 ;;
    *) die "LSF has no record of job ${jid}: it finished more than CLEAN_PERIOD ago, or the ID is wrong. Verify its outputs on disk and drop it from the graph, or pass --after-array." ;;
  esac
done

if ((${#after_ids[@]} > 0)); then
  ext_tokens=""
  for jid in "${after_ids[@]}"; do ext_tokens="${ext_tokens:+${ext_tokens},}@${jid}"; done
  for id in "${ORDER[@]}"; do
    has_internal=0
    while read -r dep; do
      [[ $dep == @* ]] || has_internal=1
    done < <(split_csv "${DEPS[$id]}")
    if ((has_internal)); then continue; fi
    if [[ ${DEPS[$id]} == "-" || -z ${DEPS[$id]} ]]; then
      DEPS[$id]=$ext_tokens
    else
      DEPS[$id]="${DEPS[$id]},${ext_tokens}"
    fi
  done
fi

# Validate everything before anything is submitted.
for id in "${ORDER[@]}"; do
  [[ -r ${SCRIPT[$id]} ]] || die "stage '${id}': cannot read '${SCRIPT[$id]}'"
  JNAME[$id]=$(jname_of "$id")
  [[ -n ${JNAME[$id]} ]] || warn "stage '${id}': no '#BSUB -J' found; treating it as a single job"
  if [[ ${JNAME[$id]} == *\[*\]* ]]; then ISARRAY[$id]=1; else ISARRAY[$id]=0; fi
  while read -r logdir; do
    [[ -d $logdir ]] ||
      warn "stage '${id}': log directory '${logdir}' does not exist; LSF will fail the job"
  done < <(awk '$1 == "#BSUB" { for (i = 2; i <= NF; i++) if ($i == "-o" || $i == "-e") print $(i + 1) }' \
    "${SCRIPT[$id]}" | xargs -r -n1 dirname | sort -u)
done

if [[ -z $state_file ]]; then
  base=${graph_file:-${spec_file:-pipeline}}
  base=${base##*/}
  base=${base%.*}
  state_file="logs/lsf-pipeline/${base}-$(date +%Y%m%d-%H%M%S).tsv"
fi

printf '\nSubmission plan (%d stages)\n\n' "${#ORDER[@]}"
for id in "${ORDER[@]}"; do
  printf '  %-16s %-44s %s\n' "$id" "${SCRIPT[$id]}" "${JNAME[$id]}"
  if [[ ${DEPS[$id]} == "-" || -z ${DEPS[$id]} ]]; then
    printf '      waits for: nothing, dispatches immediately\n'
    continue
  fi
  while read -r dep; do
    if [[ $dep == @* ]]; then
      printf '      waits for: job %s, already queued (%s)\n' \
        "${dep#@}" "$(kind_word "${EXT_ISARRAY[${dep#@}]:-0}")"
    else
      printf '      waits for: %s, every element of %s (%s)\n' \
        "$dep" "${JNAME[$dep]}" "$(kind_word "${ISARRAY[$dep]}")"
    fi
  done < <(split_csv "${DEPS[$id]}")
done

if ((dry_run)); then
  printf '\n--dry-run: nothing submitted. The state file would be %s\n' "$state_file"
  if ((want_mermaid)); then
    for id in "${ORDER[@]}"; do
      CLASS[$id]=todo
      NOTE[$id]="not submitted"
    done
    printf '\n'
    mermaid_block
  fi
  exit 0
fi

mkdir -p "$(dirname "$state_file")"
printf '#stage\tscript\tjob_name\tjob_id\tdepends_on\tdep_expr\n' >"$state_file"

printf '\nSubmitting\n\n'
for id in "${ORDER[@]}"; do
  dep_expr=""
  while read -r dep; do
    if [[ $dep == @* ]]; then
      term=$(dep_term "${dep#@}" "${EXT_ISARRAY[${dep#@}]:-0}")
    else
      [[ -n ${JOBID[$dep]:-} ]] || die "internal: stage '${dep}' has no job ID yet"
      term=$(dep_term "${JOBID[$dep]}" "${ISARRAY[$dep]}")
    fi
    dep_expr="${dep_expr:+${dep_expr} && }${term}"
  done < <(split_csv "${DEPS[$id]}")

  declare -a bsub_args=()
  if [[ -n $dep_expr ]]; then bsub_args+=(-w "$dep_expr"); fi
  declare -a extra_args=()
  read -r -a extra_args <<<"${EXTRA[$id]:-}" || true
  bsub_args+=(${extra_args[@]+"${extra_args[@]}"})

  if ! out=$(bsub ${bsub_args[@]+"${bsub_args[@]}"} <"${SCRIPT[$id]}" 2>&1); then
    printf '%s\n' "$out" >&2
    printf 'already submitted:' >&2
    for d in "${ORDER[@]}"; do
      [[ -z ${JOBID[$d]:-} ]] || printf ' %s=%s' "$d" "${JOBID[$d]}" >&2
    done
    printf '\n' >&2
    die "bsub failed for stage '${id}'; the rest of the graph was not submitted"
  fi
  jid=$(sed -n 's/^Job <\([0-9][0-9]*\)>.*/\1/p' <<<"$out" | head -n1)
  if [[ -z $jid ]]; then
    printf '%s\n' "$out" >&2
    die "could not parse a job ID out of the bsub output for stage '${id}'"
  fi
  JOBID[$id]=$jid

  printf '  %-16s %-12s %s\n' "$id" "$jid" "${dep_expr:-<no dependency>}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$id" "${SCRIPT[$id]}" "${JNAME[$id]}" "$jid" "${DEPS[$id]}" "${dep_expr:--}" >>"$state_file"
done

printf '\nState file: %s\n' "$state_file"
printf 'Check with: %s --status %s\n' "$PROG" "$state_file"
printf '        or: bjobs -A\n'

if ((want_mermaid)); then
  for id in "${ORDER[@]}"; do
    CLASS[$id]=pend
    NOTE[$id]="${JOBID[$id]} PEND"
  done
  printf '\nPaste this into the PROGRESS.md resume block:\n\n'
  mermaid_block
fi
