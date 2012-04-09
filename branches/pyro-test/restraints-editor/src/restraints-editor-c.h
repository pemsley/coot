
#ifndef RESTRAINTS_EDITOR_C_H
#define RESTRAINTS_EDITOR_C_H

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif
#endif

BEGIN_C_DECLS

void apply_restraint_by_widget(GtkWidget *w);
void restraints_editor_delete_restraint_by_widget(GtkWidget *w);
void restraints_editor_add_restraint_by_widget(GtkWidget *w);

END_C_DECLS

#endif /* RESTRAINTS_EDITOR_C_H */
