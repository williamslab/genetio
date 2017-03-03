// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include "nuclearfamily.h"

NuclearFamily::fam_ht NuclearFamily::_families;

// Given two <parents>, looks up or creates a corresponding NuclearFamily and
// adds <child> to their list of children
void NuclearFamily::addChild(PersonBulk *parents[2], PersonBulk *child) {
  par_pair_real parentsKey(parents[0], parents[1]);

  fam_ht_iter it = _families.find( &parentsKey );
  if (it == _families.end()) {
    // No family yet for these parents; create and insert it
    NuclearFamily *newFamily = new NuclearFamily();
    newFamily->_parents = new par_pair_real(parents[0], parents[1]);
    newFamily->_children.append(child); // add first child
    par_pair parentsKey = newFamily->_parents;
    _families[ parentsKey ] = newFamily;
  }
  else {
    // Have a family: append child to list of children
    // <*it> is a pair with <it->first> being the key and <it->second> the value
    // so <it->second> is the NuclearFamily
    it->second->_children.append(child);
  }
}
