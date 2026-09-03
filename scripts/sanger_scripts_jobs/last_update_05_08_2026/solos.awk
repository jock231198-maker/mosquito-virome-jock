function flush(   k, c, only) {
  if (!nm || replen < MIN) return
  c = 0; for (k in tags) { c++; only = k }
  if (c == 1) printf "%-10s %-12s %8d  %s\n", muestra, only, replen, repname
}
FNR==1 { flush(); nm=0; delete tags; replen=0
         n=split(FILENAME,p,"/"); muestra=p[n]; sub(/\.clstr$/,"",muestra) }
/^>Cluster/ { flush(); delete tags; nm=0; replen=0; repname=""; next }
{
  if (!match($0, />[A-Za-z0-9-]+__/)) next
  tags[substr($0, RSTART+1, RLENGTH-3)] = 1; nm++
  if ($0 ~ /\*[[:space:]]*$/) {
    if (match($0, /[0-9]+nt,/)) replen = substr($0, RSTART, RLENGTH-3) + 0
    if (match($0, />[^ ]+/)) { repname = substr($0, RSTART+1, RLENGTH-1); sub(/\.\.\.$/, "", repname) }
  }
}
END { flush() }
