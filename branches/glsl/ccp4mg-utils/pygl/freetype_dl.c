#ifdef _USE_DL_FREETYPE_
#include <dlfcn.h>
#include <stdio.h>

#include "freetype_dl.h"

static void *ft_lib_handle;

int init_ft(const char *filename){
  const char *error;
  error = dlerror();

  if(ft_lib_handle)
    return 1;

  ft_lib_handle = dlopen (filename, RTLD_LAZY);
  if(!ft_lib_handle){
    if ((error = dlerror()) != NULL) 
      fprintf (stderr, "%s\n", error);
    return 0;
  }
#ifdef USE_FREETYPE
  FT_Set_Pixel_Sizes= dlsym(ft_lib_handle,"FT_Set_Pixel_Sizes");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Kerning= dlsym(ft_lib_handle,"FT_Get_Kerning");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Init_FreeType= dlsym(ft_lib_handle,"FT_Init_FreeType");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_New_Face= dlsym(ft_lib_handle,"FT_New_Face");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Sfnt_Name_Count= dlsym(ft_lib_handle,"FT_Get_Sfnt_Name_Count");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Sfnt_Name= dlsym(ft_lib_handle,"FT_Get_Sfnt_Name");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Sfnt_Table= dlsym(ft_lib_handle,"FT_Get_Sfnt_Table");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Done_Face= dlsym(ft_lib_handle,"FT_Done_Face");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Postscript_Name= dlsym(ft_lib_handle,"FT_Get_Postscript_Name");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Attach_File= dlsym(ft_lib_handle,"FT_Attach_File");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Char_Index= dlsym(ft_lib_handle,"FT_Get_Char_Index");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Load_Glyph= dlsym(ft_lib_handle,"FT_Load_Glyph");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Get_Glyph= dlsym(ft_lib_handle,"FT_Get_Glyph");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Glyph_To_Bitmap= dlsym(ft_lib_handle,"FT_Glyph_To_Bitmap");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
  FT_Set_Charmap= dlsym(ft_lib_handle,"FT_Set_Charmap");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
#if defined (__APPLE_CC__)
  FT_New_Face_From_FOND= dlsym(ft_lib_handle,"FT_New_Face_From_FOND");
  if ((error = dlerror()) != NULL)  {
      fprintf (stderr, "%s\n", error);
      return 0;
  }
#endif
#endif // USE_FREETYPE
  return 1;
}
#endif // _USE_DL_
