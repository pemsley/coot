
// Triangle extras

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif


BEGIN_C_DECLS

// triangle things
int tri_count = 0; 

struct point3d {
   float x[3];
};

struct my_TRIANGLE { 
   struct point3d point[3]; 
};

#define MAX_TRIANGLES 1000000

struct my_TRIANGLE triangle_list[MAX_TRIANGLES];

void visible(int vis); 

static void idle(void); 

void inc_molecule_rot_angle(void);

END_C_DECLS

// Don't adjust this, add this information to graphics_info.
// and get rid of molecule_rot_t (or at least encapsulate 
// it in graphics_info_t.
//
class molecule_rot_t { 

 public:   
   static float x_axis_angle;
   static float y_axis_angle;
}; 
   
