#ifndef _BAM_READER_H_
#define _BAM_READER_H_

#include "reader.h"

class BAMReader : public Reader {
public:
    BAMReader(char* fname);
};

#endif /* _BAM_READER_H_ */
