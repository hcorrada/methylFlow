FILE(REMOVE_RECURSE
  "CMakeFiles/update-external-tags"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/update-external-tags.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
