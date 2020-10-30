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
#include "BKE_screen.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"

#include "bmesh.h"
#include "bmesh_tools.h"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"
#include "MOD_util.h"

static void calc_moebius_transform(float hRot[4][4], float co[3], float p)
{
  float norm;
  float h[4];
  float radius;
  int i;
  float eps;
  eps = 0.0000000000000001;

  radius = 1.0f;

  mul_v3_fl(co, 1.0f / radius);
  norm = pow(co[0], p) + pow(co[1], p) + pow(co[2], p);

  mul_v3_v3fl(h, co, 2);

  h[3] = norm - 1.0f;

  mul_v4_fl(h, 1.0f / max_ff(1.0f + norm, eps));
  mul_m4_v4(hRot, h);
  copy_v3_v3(co, h);
  mul_v3_fl(co, radius / max_ff(1.0f - h[3], eps));
}

static void moebius_transform_verts(Object *control,
                                    Object *originObject,
                                    Object *target,
                                    float (*vertexCos)[3],
                                    int numVerts,
                                    bool relocalize,
                                    float norm_power)
{
  float imat[4];
  float leftQ[4];
  float rightQ[4];
  float hRot[4][4];
  float leftMat[4][4];
  float rightMat[4][4];

  float origin[3];
  float transformedTargetPosition[3];
  int a;

  invert_m4_m4(target->imat, target->obmat);

  mat4_to_quat(leftQ, control->obmat);
  mat4_to_quat(rightQ, control->obmat);

  leftMat[0][0] = leftMat[1][1] = leftMat[2][2] = leftMat[3][3] = leftQ[0];
  leftMat[0][1] = leftMat[2][3] = -leftQ[1];
  leftMat[0][2] = leftMat[3][1] = -leftQ[2];
  leftMat[0][3] = leftMat[1][2] = -leftQ[3];
  leftMat[2][1] = leftMat[3][0] = leftQ[3];
  leftMat[1][3] = leftMat[2][0] = leftQ[2];
  leftMat[1][0] = leftMat[3][2] = leftQ[1];

  rightMat[0][0] = rightMat[1][1] = rightMat[2][2] = rightMat[3][3] = rightQ[0];
  rightMat[0][1] = rightMat[3][2] = -rightQ[1];
  rightMat[0][2] = rightMat[1][3] = -rightQ[2];
  rightMat[0][3] = rightMat[2][1] = -rightQ[3];
  rightMat[1][2] = rightMat[3][0] = rightQ[3];
  rightMat[3][1] = rightMat[2][0] = rightQ[2];
  rightMat[1][0] = rightMat[2][3] = rightQ[1];

  mul_m4_m4m4(hRot, rightMat, leftMat);

  if (relocalize) {
    zero_v3(transformedTargetPosition);
    calc_moebius_transform(hRot, transformedTargetPosition, norm_power);
  }

  if (originObject != NULL) {
    copy_v3_v3(origin, originObject->obmat[3]);
  }
  else {
    copy_v3_v3(origin, control->obmat[3]);
  }
  for (a = 0; a < numVerts; a++) {
    mul_m4_v3(target->obmat, vertexCos[a]);
    sub_v3_v3(vertexCos[a], origin);

    calc_moebius_transform(hRot, vertexCos[a], norm_power);

    add_v3_v3(vertexCos[a], origin);
    mul_m4_v3(target->imat, vertexCos[a]);

    if (relocalize) {
      sub_v3_v3(vertexCos[a], transformedTargetPosition);
    }
  }
}

static void deformVerts(ModifierData *md,
                        const ModifierEvalContext *ctx,
                        Mesh *mesh,
                        float (*vertexCos)[3],
                        int numVerts)
{
  MoebiusModifierData *mmd = (MoebiusModifierData *)md;
  moebius_transform_verts(mmd->control,
                          mmd->origin,
                          ctx->object,
                          vertexCos,
                          numVerts,
                          mmd->flags & eMoebiusModifierFlag_localize,
                          mmd->norm_power);
}

static void deformVertsEM(ModifierData *md,
                          const ModifierEvalContext *ctx,
                          struct BMEditMesh *editData,
                          Mesh *mesh,
                          float (*vertexCos)[3],
                          int numVerts)
{
  MoebiusModifierData *mmd = (MoebiusModifierData *)md;
  moebius_transform_verts(mmd->control,
                          mmd->origin,
                          ctx->object,
                          vertexCos,
                          numVerts,
                          mmd->flags & eMoebiusModifierFlag_localize,
                          mmd->norm_power);
}

/* Moebius Transform */
static void initData(ModifierData *md)
{
  MoebiusModifierData *smd = (MoebiusModifierData *)md;
  smd->norm_power = 2.0;
}

static void copyData(ModifierData *md, ModifierData *target)
{
  MoebiusModifierData *smd = (MoebiusModifierData *)md;
  MoebiusModifierData *tsmd = (MoebiusModifierData *)target;

  tsmd->control = smd->control;
  tsmd->origin = smd->origin;
  tsmd->flags = smd->flags;
  tsmd->norm_power = smd->norm_power;
}

static void foreachIDLink(ModifierData *md, Object *ob, IDWalkFunc walk, void *userData)
{
  MoebiusModifierData *mmd = (MoebiusModifierData *)md;

  walk(userData, ob, (ID **)&mmd->control, IDWALK_NOP);
  walk(userData, ob, (ID **)&mmd->origin, IDWALK_NOP);
}

static bool isDisabled(const Scene *UNUSED(scene), ModifierData *md, bool UNUSED(userRenderParams))
{
  MoebiusModifierData *mmd = (MoebiusModifierData *)md;

  return !mmd->control;
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
  MoebiusModifierData *lmd = (MoebiusModifierData *)md;

  if (lmd->control)
    DEG_add_object_relation(
        ctx->node, lmd->control, DEG_OB_COMP_TRANSFORM, "Moebius Transformation Modifier");
  if (lmd->origin)
    DEG_add_object_relation(
        ctx->node, lmd->origin, DEG_OB_COMP_TRANSFORM, "Moebius Transformation Modifier");
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA ob_ptr;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, &ob_ptr);

  uiLayoutSetPropSep(layout, true);

  uiItemR(layout, ptr, "control", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "origin", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "localize", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "norm_power", 0, NULL, ICON_NONE);

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  modifier_panel_register(region_type, eModifierType_Moebius, panel_draw);
}

ModifierTypeInfo modifierType_Moebius = {
    /* name */ "Moebius Transformation",
    /* structName */ "MoebiusModifierData",
    /* structSize */ sizeof(MoebiusModifierData),
    /* srna */ &RNA_MoebiusModifier,
    /* type */ eModifierTypeType_Constructive,
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
