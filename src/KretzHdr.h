#ifndef KretzHdr_h
#define KretzHdr_h

namespace itk

struct Tag
{
	unsigned short group;
	unsigned short element;
};

std::map<std::string, Tag> tag_dict;
tag_dict.insert(std::make_pair("PatientName", Tag(0x0110, 0x0002) ));
tag_dict.insert(std::make_pair("Dimension_x", Tag(0xc000, 0x0001) ));

enum Kretz 
{
  Tag Patient_Name = {0x0110, 0x0002}
  Tag Dimension_x = {0xc000, 0x0001}
  Tag Dimension_y = {0xc000, 0x0002}
  Tag Dimension_z = {0xc000, 0x0003}
  Tag Resolution = {0xc100, 0x0001}
  Tag Offset1 = {0xc200, 0x0001}
  Tag Offset2 = {0xc200, 0x0002}
  Tag Angles1 = {0xc300, 0x0001}
  Tag Angles2 = {0xc300, 0x0002}
  Tag ImageData = {0xd000, 0x0001}
  Tag CineFrames = {0xd400, 0x0001}
  Tag SizeOfFrame = {0xd400, 0x0002}
  Tag TimingOfFrame = {0xd400, 0x0005}
  Tag ImageData4D = {0xd600, 0x0001}
};


#endif  /* __KretzHdr_h */
