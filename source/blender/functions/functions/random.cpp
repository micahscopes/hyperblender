#include "FN_types.hpp"
#include "FN_tuple_call.hpp"

#include "BLI_lazy_init_cxx.h"
#include "BLI_math_bits.h"

#include "random.hpp"

namespace FN {
namespace Functions {

using namespace Types;

static uint32_t random_int(uint32_t x)
{
  x = (x << 13) ^ x;
  return x * (x * x * 15731 + 789221) + 1376312589;
}

static float random_float(uint32_t x)
{
  x = random_int(x);
  return (float)x / 4294967296.0f;
}

class RandomNumber : public TupleCallBody {
  void call(Tuple &fn_in, Tuple &fn_out, ExecutionContext &UNUSED(ctx)) const override
  {
    FN_TUPLE_CALL_NAMED_REF(this, fn_in, fn_out, inputs, outputs);

    float seed = inputs.get<float>(0, "Seed");
    float min = inputs.get<float>(1, "Min");
    float max = inputs.get<float>(2, "Max");
    float result = random_float(float_as_uint(seed)) * (max - min) + min;
    outputs.set<float>(0, "Value", result);
  }
};

BLI_LAZY_INIT(SharedFunction, GET_FN_random_number)
{
  FunctionBuilder builder;
  builder.add_input("Seed", TYPE_float);
  builder.add_input("Min", TYPE_float);
  builder.add_input("Max", TYPE_float);
  builder.add_output("Value", TYPE_float);

  auto fn = builder.build("Random Number");
  fn->add_body<RandomNumber>();
  return fn;
}

}  // namespace Functions
}  // namespace FN
