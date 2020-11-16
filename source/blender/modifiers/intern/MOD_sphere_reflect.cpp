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

static void bm_sphere_reflect(Object *sphereObject, Object *target, BMesh *bm)
{
  auto sphere = sphere_of_object(sphereObject);
  int j;
  double refCos[3];

  BMVert *v;
  BMIter iter;
  BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
    mul_m4_v3(target->obmat, v->co);
    for (j = 0; j < 3; j++)
      refCos[j] = (double)(v->co)[j];
    auto p = point(refCos[0], refCos[1], refCos[2]);
    auto ref = apply_even_versor(sphere, p);
    auto v_ref = -ref / abs(ref | ni);
    refCos[0] = v_ref | e1;
    refCos[1] = v_ref | e2;
    refCos[2] = v_ref | e3;
    for (j = 0; j < 3; j++)
      (v->co)[j] = (float)refCos[j];
    mul_m4_v3(target->imat, v->co);
  }
}

static Mesh *adaptiveRefine(SphereReflectModifierData *smd,
                            const ModifierEvalContext *ctx,
                            Mesh *mesh)
{
  Mesh *result;
  BMesh *bm;
  const struct BMeshCreateParams bm_create_params = {
      .use_toolflags = true,
  };
  const struct BMeshFromMeshParams bm_from_mesh_params = {0};
  bm = BKE_mesh_to_bmesh_ex(mesh, &bm_create_params, &bm_from_mesh_params);

  if (smd->flags & MOD_SPHERE_REFLECT_ADAPTIVE) {
    BMOperator op;
    BMOpSlot *slot_geom_inner;
    BMElem **subdivided_elem_array;
    BMVert *v;
    BMFace *f;
    BMEdge *e;
    BMIter iter;
    BMIter fiter;
    BMIter eiter;
    float max_distance = 1.0;
    int num_close_verts = -1;
    int cuts = smd->cuts;
    int i;
    float worldCo[3];
    for (i = 0; i < smd->iterations; i++) {
      // printf("iteration %i\n", i);
      BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
      max_distance = 1.0 * pow(1.0 / (float(cuts) + 1.0), i) / smd->falloff;

      if (i == 0) {
        num_close_verts = 0;
        BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
          copy_v3_v3(worldCo, v->co);
          mul_m4_v3(ctx->object->obmat, worldCo);
          if (len_v3v3(worldCo, smd->sphere->obmat[3]) < max_distance) {
            num_close_verts += 1;
            BM_ITER_ELEM (f, &fiter, v, BM_FACES_OF_VERT) {
              BM_ITER_ELEM (e, &eiter, f, BM_EDGES_OF_FACE) {
                BM_elem_flag_enable(e, BM_ELEM_TAG);
              }
            }
          }
        }
      }
      else if (num_close_verts == 0) {
        // printf("num_close_verts == 0\n");
        break;
      }
      else {
        num_close_verts = 0;
        slot_geom_inner = BMO_slot_get(op.slots_out, "geom_inner.out");
        subdivided_elem_array = (BMElem **)slot_geom_inner->data.buf;
        // printf("amount last subdivided %i\n", slot_geom_inner->len);
        for (int j = 0; j < slot_geom_inner->len; j++) {
          if (subdivided_elem_array[i]->head.htype == BM_VERT) {
            // printf("checking BM_VERT %i\n", j);
            v = (BMVert *)subdivided_elem_array[j];
            copy_v3_v3(worldCo, v->co);
            mul_m4_v3(ctx->object->obmat, worldCo);
            if (len_v3v3(worldCo, smd->sphere->obmat[3]) < max_distance) {
              num_close_verts += 1;
              BM_ITER_ELEM (f, &fiter, v, BM_FACES_OF_VERT) {
                BM_ITER_ELEM (e, &eiter, f, BM_EDGES_OF_FACE) {
                  BM_elem_flag_enable(e, BM_ELEM_TAG);
                }
              }
            }
          }
        }
      }
      BMO_op_initf(bm,
                   &op,
                   BMO_FLAG_DEFAULTS,
                   "subdivide_edges edges=%he cuts=%i use_grid_fill=%b",
                   BM_ELEM_TAG,
                   cuts,
                   true);
      BMO_op_exec(bm, &op);
    }
  }

  bm_sphere_reflect(smd->sphere, ctx->object, bm);
  result = BKE_mesh_from_bmesh_for_eval_nomain(bm, NULL, mesh);
  BM_mesh_free(bm);

  return result;
}

static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *ctx, Mesh *mesh)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
  Mesh *result;
  if (!(result = adaptiveRefine(smd, ctx, mesh))) {
    return mesh;
  }

  return result;
}

// static void deformVerts(ModifierData *md,
//                         const ModifierEvalContext *ctx,
//                         Mesh *mesh,
//                         float (*vertexCos)[3],
//                         int numVerts)
// {
//   SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
//   sphere_reflect(smd->sphere, ctx->object, vertexCos, numVerts);
// }

// static void deformVertsEM(ModifierData *md,
//                           const ModifierEvalContext *ctx,
//                           struct BMEditMesh *editData,
//                           Mesh *mesh,
//                           float (*vertexCos)[3],
//                           int numVerts)
// {
//   SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
//   sphere_reflect(smd->sphere, ctx->object, vertexCos, numVerts);
// }

/* SphereReflect Transform */
static void initData(ModifierData *md)
{
  SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
  smd->theta = 0.0;
  smd->iterations = 8;
  smd->cuts = 1;
  smd->falloff = 2.0;
}

// static void copyData(ModifierData *md, ModifierData *target)
// {
//   SphereReflectModifierData *smd = (SphereReflectModifierData *)md;
//   SphereReflectModifierData *tsmd = (SphereReflectModifierData *)target;

//   tsmd->sphere = smd->sphere;
// }

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
  uiItemR(layout, ptr, "adaptive", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "cuts", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "iterations", 0, NULL, ICON_NONE);
  uiItemR(layout, ptr, "falloff", 0, NULL, ICON_NONE);

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
    /* type */ eModifierTypeType_Constructive,
    /* flags */
    ModifierTypeFlag(eModifierTypeFlag_AcceptsMesh | eModifierTypeFlag_SupportsEditmode |
                     eModifierTypeFlag_AcceptsCVs),
    /* icon */ NULL, /* TODO: Use correct icon. */
    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ NULL,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ NULL,
    /* deformMatricesEM */ NULL,
    /* modifyMesh */ modifyMesh,
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
