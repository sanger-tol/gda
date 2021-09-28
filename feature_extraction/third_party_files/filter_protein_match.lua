name = "Protein Domain Filter"
author = "Sascha Kastens"
version = "1.0"
email = "mail@skastens.de"
short_descr = "Filters out candidates without protein domains"
description = "Filters out a candidate if it does not contain at " ..
              "least one node of type 'protein_match'."

function filter(gn)
  gfi = gt.feature_node_iterator_new(gn)
  node = gfi:next()
  while not (node == nil) do
    if (node:get_type() == "protein_match") then
      return false
    end
    node = gfi:next()
  end
  return true
end
