// Library for I/O of genetic data
// Author: Amy Williams <alw289  cornell edu>
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <stdio.h>
#include "superperson.h"


////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<const char *> SuperPerson::_popLabels;
int SuperPerson::_maxPersonIdLength = 0;
int SuperPerson::_numMales = 0;
int SuperPerson::_numSexUnknown = 0;

SuperPerson::SuperPerson(char *id, char sex, int popIndex,
			 short familyIdLength) {
  _popIndex = popIndex;
  _ignore = (_popIndex == -1) ? true : false;
  if (!_ignore) {
    int idLength = strlen(id);
    if (idLength > MAX_PERSON_ID - 1) { // -1 for null terminating character
      fprintf(stderr, "Error: id %s is more than the maximum length of %d characters\n",
	      id, MAX_PERSON_ID - 1);
      exit(5);
    }
    if (idLength > _maxPersonIdLength)
      _maxPersonIdLength = idLength;
    strcpy(_id, id);
    _familyIdLength = familyIdLength;
    _sex = sex;

    if (_sex == 'M')
      _numMales++;
    else if (_sex == 'U')
      _numSexUnknown++;
    else if (_sex != 'F') {
      fprintf(stderr, "\nERROR: unknown sex %c for id %s!\n", _sex,
	      _id);
      exit(4);
    }

  }
}

SuperPerson::~SuperPerson() {
}
