#include "GAMER.h"
#include "NuclearEos.h"

void Init_NuclearEos()
{
  nuc_eos_C_ReadTable("LS220.h5");   // TODO: set table name as a runtime parameter

  if (MPI_Rank == 0)
  {
    printf("nrho: %d\n",nrho);
    printf("ntemp: %d\n",ntemp);
    printf("nye: %d\n",nye);
  }

  //nuc_eos_C_testing();

}
