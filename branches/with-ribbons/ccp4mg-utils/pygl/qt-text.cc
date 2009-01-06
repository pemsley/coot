/* basically dummy file as we dont want to use qt! */
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include "ppmutil.h"
#include "cprimitive.h"

/*
#ifdef __WIN32__
#include <QApplication>
static QApplication *app = 0;
#endif
*/

void SimpleText::renderStringToPixmap(){

/*
#ifdef __WIN32__
  // Why is this necessary?
  if(app==0){
    int argc = 0;
    char **argv = 0;
    app = new QApplication(argc,argv);
  }
#endif
*/
  std::vector<std::string> strings(1);
  strings[0] = text; // For the moment.
  
  //QFont font(family.c_str(),fn_size); // Need to consider slant/weight
  //if(slant=="i") font.setItalic(true);
  //if(slant=="o") font.setItalic(true);
  //if(weight=="bold") font.setBold(true);
  //QFontMetrics fm(font);

  int max_width = 0;
  //for(unsigned i=0;i<strings.size();i++){
  //  int width = fm.boundingRect(strings[i].c_str()).width();
  //  if(width>max_width) max_width = width;
  //}

  int lines = strings.size();

  int required_width = max_width;
  //int required_height = fm.leading()*(lines-1)+fm.height()*lines;

  //required_width = 2<<int(ceil(log(required_width)/log(2))-1);
  //required_height =  2<<int(ceil(log(required_height)/log(2))-1);

  //QImage *pixmap = new QImage(required_width,required_height,QImage::Format_ARGB32);
  //pixmap->fill(0);

  //QPainter painter;
  //painter.begin(pixmap);
  //painter.setFont(font);
  //painter.setPen(QColor(255,255,255,255));
  //for(unsigned i=0;i<strings.size();i++){
  //  painter.drawText(-fm.boundingRect(strings[i].c_str()).left(),fm.height()-fm.descent()+i*fm.lineSpacing(),strings[i].c_str());
 // }
  //painter.end();

  //texture = image_info(pixmap->width(),pixmap->height(),pixmap->bits(),IMAGEINFO_RGBA);
  //texture.invert();

  // Make this new!
  texture_id = 0;

  //delete pixmap;

}
