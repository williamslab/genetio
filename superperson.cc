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
int SuperPerson::_numGenderUnknown = 0;

SuperPerson::SuperPerson(char *id, char gender, int popIndex,
			 short familyIdLength) {
  _popIndex = popIndex;
  _ignore = (_popIndex == -1) ? true : false;
  if (!_ignore) {
    int idLength = strlen(id);
    if (idLength > _maxPersonIdLength)
      _maxPersonIdLength = idLength;
    _id = new char[ idLength + 1 ];
    strcpy(_id, id);
    _familyIdLength = familyIdLength;
    _gender = gender;

    if (_gender == 'M')
      _numMales++;
    else if (_gender == 'U')
      _numGenderUnknown++;
    else if (_gender != 'F') {
      fprintf(stderr, "\nERROR: unknown gender %c for id %s!\n", _gender,
	      _id);
      exit(4);
    }

  }
}

SuperPerson::~SuperPerson() {
  if (!_ignore)
    delete [] _id;
}
