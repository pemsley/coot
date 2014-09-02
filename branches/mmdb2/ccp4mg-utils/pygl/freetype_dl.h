#ifndef _FREETYPEDL_H_
#define _FREETYPEDL_H_

#ifdef USE_FREETYPE

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_FREETYPE_H
#include FT_ERRORS_H
#include FT_MULTIPLE_MASTERS_H                                   
#include FT_GLYPH_H                                              
#include FT_SFNT_NAMES_H
#include FT_TRUETYPE_IDS_H
#include FT_TRUETYPE_TABLES_H

#if defined (__APPLE_CC__)
#include <Carbon/Carbon.h>
#include FT_MAC_H
#endif

#ifdef __cplusplus
extern "C"{
#endif
int init_ft(const char *filename);

FT_Error (*FT_Set_Pixel_Sizes_fun)( FT_Face face, FT_UInt pixel_width, FT_UInt pixel_height );
FT_Error (*FT_Get_Kerning_fun)( FT_Face face, FT_UInt left_glyph, FT_UInt right_glyph, FT_UInt kern_mode, FT_Vector *akerning );
FT_Error (*FT_Init_FreeType_fun)( FT_Library  *alibrary );
FT_Error (*FT_New_Face_fun)( FT_Library   library, const char* filepathname, FT_Long face_index, FT_Face *aface ); 
FT_UInt (*FT_Get_Sfnt_Name_Count_fun)( FT_Face  face );
FT_Error (*FT_Get_Sfnt_Name_fun)( FT_Face face, FT_UInt idx, FT_SfntName *aname );
void* (*FT_Get_Sfnt_Table_fun)( FT_Face face, FT_Sfnt_Tag tag );
FT_Error (*FT_Done_Face_fun)( FT_Face face );
const char* (*FT_Get_Postscript_Name_fun)( FT_Face  face );
FT_Error (*FT_Attach_File_fun)( FT_Face face, const char* filepathname );
FT_UInt (*FT_Get_Char_Index_fun)( FT_Face face, FT_ULong charcode );
FT_Error (*FT_Load_Glyph_fun)( FT_Face face, FT_UInt glyph_index, FT_Int32 load_flags );
FT_Error (*FT_Get_Glyph_fun)( FT_GlyphSlot slot, FT_Glyph *aglyph );
FT_Error (*FT_Glyph_To_Bitmap_fun)( FT_Glyph* the_glyph, FT_Render_Mode render_mode, FT_Vector* origin, FT_Bool destroy );
FT_Error (*FT_Set_Charmap_fun)(FT_Face face, FT_CharMap charmap);
#if defined (__APPLE_CC__)
FT_Error (*FT_New_Face_From_FOND_fun)(FT_Library library, Handle fond, FT_Long face_index, FT_Face *aface );
#endif
#ifdef __cplusplus
}
#endif

#define FT_Set_Pixel_Sizes FT_Set_Pixel_Sizes_fun
#define FT_Get_Kerning FT_Get_Kerning_fun
#define FT_Init_FreeType FT_Init_FreeType_fun
#define FT_New_Face FT_New_Face_fun
#define FT_Get_Sfnt_Name_Count FT_Get_Sfnt_Name_Count_fun
#define FT_Get_Sfnt_Name FT_Get_Sfnt_Name_fun
#define FT_Get_Sfnt_Table FT_Get_Sfnt_Table_fun
#define FT_Done_Face FT_Done_Face_fun
#define FT_Get_Postscript_Name FT_Get_Postscript_Name_fun
#define FT_Attach_File FT_Attach_File_fun
#define FT_Get_Char_Index FT_Get_Char_Index_fun
#define FT_Load_Glyph FT_Load_Glyph_fun
#define FT_Get_Glyph FT_Get_Glyph_fun
#define FT_Glyph_To_Bitmap FT_Glyph_To_Bitmap_fun
#define FT_Set_Charmap FT_Set_Charmap_fun
#if defined (__APPLE_CC__)
#define FT_New_Face_From_FOND FT_New_Face_From_FOND_fun
#endif

#endif
#endif // _FREETYPEDL_H_
