
#include <names_map.h>

// Add a name and associate it with an index
void NamesMap::add(const std::string &name, int index) {
  name_to_index[name] = index;
  index_to_name[index] = name;
}

// Get the index associated with a name
int NamesMap::get_index(const std::string &name) const {
  auto it = name_to_index.find(name);
  return (it != name_to_index.end()) ? it->second : -1; // -1 for not found
}

// Get the name associated with an index
std::string NamesMap::get_name(int index) const {
  auto it = index_to_name.find(index);
  return (it != index_to_name.end()) ? it->second
                                     : ""; // Empty string for not found
}

// Get a vector of all names
std::vector<std::string> NamesMap::get_all_names() const {
  std::vector<std::string> names;
  for (const auto &pair : index_to_name) {
    names.push_back(pair.second);
  }
  return names;
}
