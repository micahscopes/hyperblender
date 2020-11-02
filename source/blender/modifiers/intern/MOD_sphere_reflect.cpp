#include <string.h>

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLT_translation.h"

#include "DNA_defaults.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_screen_types.h"

#include "BKE_context.h"
#include "BKE_lib_id.h"
#include "BKE_lib_query.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_screen.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"

#include "bmesh.h"
#include "bmesh_tools.h"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"
extern "C" {
#include "MOD_util.h"
}

#include "gatl/ga3c.hpp"
using namespace ga3c;

auto sphere_of_object(Object *object)
{
  float o[3];
  copy_v3_v3(o, object->obmat[3]);
  float size[3] = {1.0f, 1.0f, 1.0f};
  float radius = 1.0f;
  BKE_object_dimensions_get(object, size);
  radius = max_fff(size[0], size[1], size[2]) * 0.5f;
  auto p1 = point(o[0] + radius, o[1], o[2]);
  auto p2 = point(o[0] - radius, o[1], o[2]);
  auto p3 = point(o[0], o[1] + radius, o[2]);
  auto p4 = point(o[0], o[1], o[2] + radius);
  auto sphere = p1 ^ p2 ^ p3 ^ p4;
  return sphere;
}

static void sphere_reflect(Object *sphereObject,
                           Object *target,
                           float (*vertexCos)[3],
                           int numVerts)
{
  auto sphere = sphere_of_object(sphereObject);
  const auto Mnk = ni ^ no;
  int i;
  int j;
  double refCos[3];
  for (i = 0; i < numVerts; i++) {
    mul_m4_v3(target->obmat, vertexCos[i]);
    for (j = 0; j < 3; j++)
      refCos[j] = (double)vertexCos[i][j];
    auto p = point(refCos[0], refCos[1], refCos[2]);
    auto ref = apply_even_versor(sphere, p);
    auto v_ref = -ref / abs(ref | ni);
    refCos[0] = v_ref | e1;
    refCos[1] = v_ref | e2;
    refCos[2] = v_ref | e3;
    for (j = 0; j < 3; j++)
      vertexCos[i][j] = (float)refCos[j];
    // copy_v3_v3(vertexCos[i], refCos);
    mul_m4_v3(target->imat, vertexCos[i]);
  }
}

static void deformVerts(ModifierData *md,
                        const ModifierEvalContext *ctx,
                        Mesh *mesh,
                        float (*vertexCos)[3],
                        int numVerts)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
  sphere_reflect(smd->sphere, ctx->object, vertexCos, numVerts);
}

static void deformVertsEM(ModifierData *md,
                          const ModifierEvalContext *ctx,
                          struct BMEditMesh *editData,
                          Mesh *mesh,
                          float (*vertexCos)[3],
                          int numVerts)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
  sphere_reflect(smd->sphere, ctx->object, vertexCos, numVerts);
}

/* SphereReflect Transform */
static void initData(ModifierData *md)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
}

static void copyData(ModifierData *md, ModifierData *target)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
  SphereReflectModifierData *tsmd = (SphereReflectModifierData *)target;

  tsmd->sphere = smd->sphere;
}

static void foreachIDLink(ModifierData *md, Object *ob, IDWalkFunc walk, void *userData)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;

  walk(userData, ob, (ID **)&smd->sphere, IDWALK_NOP);
}

static bool isDisabled(const Scene *UNUSED(scene), ModifierData *md, bool UNUSED(userRenderParams))
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;

  return smd->sphere == NULL;
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
  SphereReflectModifierData *lmd = (SphereReflectModifierData *)md;

  if (lmd->sphere)
    DEG_add_object_relation(
        ctx->node, lmd->sphere, DEG_OB_COMP_TRANSFORM, "SphereReflect Transformation Modifier");
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA ob_ptr;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, &ob_ptr);

  uiLayoutSetPropSep(layout, true);

  uiItemR(layout, ptr, "sphere", 0, NULL, ICON_NONE);

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  modifier_panel_register(region_type, eModifierType_SphereReflect, panel_draw);
}

ModifierTypeInfo modifierType_SphereReflect = {
    /* name */ "SphereReflect Transformation",
    /* structName */ "SphereReflectModifierData",
    /* structSize */ sizeof(SphereReflectModifierData),
    /* srna */ &RNA_SphereReflectModifier,
    /* type */ eModifierTypeType_OnlyDeform,
    /* flags */
    ModifierTypeFlag(eModifierTypeFlag_AcceptsMesh | eModifierTypeFlag_SupportsEditmode |
                     eModifierTypeFlag_AcceptsCVs),
    /* icon */ NULL, /* TODO: Use correct icon. */
    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ deformVerts,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ deformVertsEM,
    /* deformMatricesEM */ NULL,
    /* modifyMesh */ NULL,
    /* modifyHair */ NULL,
    /* modifyPointCloud */ NULL,
    /* modifyVolume */ NULL,

    /* initData */ initData,
    /* requiredDataMask */ NULL,
    /* freeData */ NULL,
    /* isDisabled */ isDisabled,
    /* updateDepsgraph */ updateDepsgraph,
    /* dependsOnTime */ NULL,
    /* dependsOnNormals */ NULL,
    /* foreachIDLink */ foreachIDLink,
    /* foreachTexLink */ NULL,
    /* freeRuntimeData */ NULL,
    /* panelRegister */ panelRegister,
    /* blendWrite */ NULL,
    /* blendRead */ NULL,
};
