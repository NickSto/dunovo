
substr($0, 1, 1) == ">" {
  header = $0
  split(header, fields1, ".")
  split(fields1[2], fields2)
  mate = fields2[1]
  if (mate == target) {
    print fields1[1]" "fields2[2]
  }
  next
}
{
  if (mate == target) {
    print
  }
}
