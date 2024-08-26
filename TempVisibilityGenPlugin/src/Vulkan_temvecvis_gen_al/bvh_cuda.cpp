
#include "g_common.h"

#include "bvh_cuda.h"

int bvh_size( my_node* root, int addr )
{
  int blk_size = int(
    sizeof(root->left) +
    sizeof(root->right) +
    sizeof(root->box.m) +
    sizeof(root->box.n) +
    sizeof(root->conec_up) +
    sizeof(root->cos_coner_up) +
    sizeof(root->conec_down) +
    sizeof(root->cos_coner_down) +
    sizeof(root->idx.size()) +
    root->idx.size() * sizeof(int));
  int curr_addr = addr+blk_size;
  if( root->left )
    curr_addr = bvh_size(root->left, curr_addr);
  if( root->right )
    curr_addr = bvh_size(root->right, curr_addr);
  return curr_addr;
}

int bvh_memcpy( my_node* root, unsigned char *bvh0, int addr )
{
  int blk_size = int(
    sizeof(root->left) +
    sizeof(root->right) +
    sizeof(root->box.m) +
    sizeof(root->box.n) +
    sizeof(root->conec_up) +
    sizeof(root->cos_coner_up) +
    sizeof(root->conec_down) +
    sizeof(root->cos_coner_down) +
    sizeof(root->idx.size()) +
    root->idx.size() * sizeof(int));

  int i;
  cubvh *t = (cubvh*) (bvh0+addr);
    t->left = 0;
    t->right = 0;
    t->m = make_float3(root->box.m.x, root->box.m.y, root->box.m.z);
    t->n = make_float3(root->box.n.x, root->box.n.y, root->box.n.z);
    t->conec_up = make_float3(root->conec_up.x, root->conec_up.y, root->conec_up.z);
    t->cos_coner_up = root->cos_coner_up;
    t->conec_down = make_float3(root->conec_down.x, root->conec_down.y, root->conec_down.z);
    t->cos_coner_down = root->cos_coner_down;
    t->ni = int(root->idx.size());
    int *idx = t->idx;//(int*)( ((unsigned char*)t)+sizeof(cubvh) );
    //printf("root->idx.size() %i\n", root->idx.size() );
    for(i=0; i<root->idx.size(); i++ )
    {
      idx[i] = root->idx[i];
      //printf("  t->idx[%i] %i\n", i, t->idx[i] );
    }
  int curr_addr = addr+blk_size;
  if( root->left )
  {
    t->left = curr_addr;
    curr_addr = bvh_memcpy(root->left, bvh0, curr_addr);
  }
  if( root->right )
   {
    t->right = curr_addr;
    curr_addr = bvh_memcpy(root->right, bvh0, curr_addr);
  }

  return curr_addr;
}

void bvh_print( unsigned char *bvh0, int addr, int &m )
{
  cubvh *ccc = (cubvh*) (bvh0+addr);
  int *idx = ccc->idx;//(int*)( ((unsigned char*)ccc)+sizeof(cubvh) );

  int i;
  printf("ccc->ni %i\n", ccc->ni );
  for(i=0; i<ccc->ni; i++ )
    printf("  ccc->idx[%i] %i\n", i, idx[i] );

  //printf( "%f %f %f\n", ccc->m.x, ccc->m.y, ccc->m.z );
  //printf( "%f %f %f\n", ccc->n.x, ccc->n.y, ccc->n.z );
  //printf( "%i\n", ccc->left );
  //printf( "%i\n", ccc->right );
  //printf( "%i\n\n", ccc->ni );

  m = m + ccc->ni;


  if( ccc->left )
    bvh_print(bvh0, ccc->left, m);
  if( ccc->right )
    bvh_print(bvh0, ccc->right, m);
}

