/**********************************************************************
The MIT License

Copyright (c) 2014 Intel Corporation

	Permission is hereby granted, free of charge, to any person
	obtaining a copy of this software and associated documentation
	files (the "Software"), to deal in the Software without
	restriction, including without limitation the rights to use,
	copy, modify, merge, publish, distribute, sublicense, and/or
	sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following
	conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the
	Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
	KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
	WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
	COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
	OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**********************************************************************/
#include "options.h"
#if ((MAJOR_VERSION == IGZIP0) || (MAJOR_VERSION == IGZIP1))

#include "types.h"

//extern "C" UINT32 CrcTable[256];

//void
//init_crc(UINT32 &crc)
//{
//    crc = -1;
//}

//void
//update_crc(UINT32 &crc, UINT8 byte)
//{
//    crc = (crc >> 8) ^ CrcTable[(crc ^ byte) & 0xFF];
//}

UINT32
get_crc(UINT32 crc)
{
    return ~crc;
}

#endif // if ((MAJOR_VERSION == IGZIP0) || (MAJOR_VERSION == IGZIP1))
