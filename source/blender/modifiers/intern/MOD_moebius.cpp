#include "BKE_modifier.h"
#include "DNA_mesh_types.h"
#include "gal/pga.hpp"

using namespace gal::pga;
using gal::compute;

static Mesh *moebius_applyModifier(struct ModifierData *md,
                                   const struct ModifierEvalContext *ctx,
                                   struct Mesh *mesh)
{
  point<> p1{1, 0, 0}; // The default value type is `float`
  point<> p2{0, 1, 1};

  // Take the regressive product of the two points to construct
  // the line passing through them.
  line<> l = compute([](auto p1, auto p2) {
          return p1 & p2;
      }, p1, p2);
  return mesh;
}

ModifierTypeInfo modifierType_Moebius = {
    /* name */ "Moebius Transformation",
    /* structName */ "MoebiusModifierData",
    /* structSize */ sizeof(MoebiusModifierData),
    /* type */ eModifierTypeType_DeformOrConstruct,
    /* flags */ 
    ModifierTypeFlag(
      eModifierTypeFlag_AcceptsMesh
      | eModifierTypeFlag_SupportsEditmode
      | eModifierTypeFlag_AcceptsCVs
    ) ,
    /* copyData */ NULL,

    /* deformVerts */ NULL,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ NULL,
    /* deformMatricesEM */ NULL,
    /* applyModifier */ moebius_applyModifier,

    /* initData */ NULL,
    /* requiredDataMask */ NULL,
    /* freeData */ NULL,
    /* isDisabled */ NULL,
    /* updateDepsgraph */ NULL,
    /* dependsOnTime */ NULL,
    /* dependsOnNormals */ NULL,
    /* foreachObjectLink */ NULL,
    /* foreachIDLink */ NULL,
    /* foreachTexLink */ NULL,
    /* freeRuntimeData */ NULL,
};
