#include <direct.h>
#include <windows.h>

#include "g_filedialog.h"
#include "g_common.h"



GFileDialog::GFileDialog()
{
  strcpy( spath, "" );
  strcpy( spath_relative, "" );
  strcpy( init_dir, "." );
}

bool GFileDialog::prompt_open( const char *file_filter )
{
  replace_char( init_dir, '/', '\\' );

  OPENFILENAMEA ofn = { 0 };
    ofn.lStructSize  = sizeof(OPENFILENAMEA);
    ofn.lpstrFilter  = file_filter;
    ofn.nFilterIndex = 1;
    ofn.lpstrFile    = spath;
    ofn.lpstrInitialDir = init_dir;
    ofn.nMaxFile     = G_FILEDIALOG_NMAXFILE;
    ofn.Flags        = OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;

  char curr_dir[256], dir0[256], dir1[256];

  _getcwd( curr_dir, 256 );
  _chdir(init_dir);
  _getcwd( dir0, 256 );
  _chdir(curr_dir);


  if( GetOpenFileNameA(&ofn) )
  {
    _getcwd( dir1, 256 );

    int a, b, c;
    int i, i0, l;
    if(dir0[strlen(dir0)-1]!='\\')
      strcat(dir0,"\\");
    if(dir1[strlen(dir1)-1]!='\\')
      strcat(dir1,"\\");
    a = strlen(dir0);
    b = strlen(dir1);
    c = g_min(a,b);
    for( i=0; i<c && dir0[i]==dir1[i]; i++ )
    {
      if( dir0[i]=='\\' )
        i0 = i;
    }
    if(i0==0)
    {
      strcpy( spath_relative, spath );
    }else
    {
      for( i=i0, l=0; i<a; i++ )
        if( dir0[i]=='\\' )
          l++;
      GPath gp = parse_spath( spath );
      strcpy( spath_relative, "" );
      for( i=1; i<l; i++ )
      strcat( spath_relative, "..\\" );
      strcat( spath_relative, &dir1[i0+1] );
      strcat( spath_relative, gp.fname );
      if(gp.ename[0])
      {
        strcat( spath_relative, "." );
        strcat( spath_relative, gp.ename );
      }
    }
    // printf( "%s\n", spath_relative );


    _chdir(curr_dir);
    return true;
  }

  return false;
}

bool GFileDialog::prompt_save( const char *file_filter )
{
  char curr_dir[256];
  _getcwd( curr_dir, 256 );
  GPath gp = parse_spath(spath);
  sprintf( spath, "%s%s", gp.dname, gp.fname );


  OPENFILENAMEA ofn ={0};
  ofn.lStructSize  = sizeof(OPENFILENAMEA);
  ofn.lpstrFilter  = file_filter;
  ofn.nFilterIndex = 1;
  ofn.nMaxFile     = G_FILEDIALOG_NMAXFILE;
  ofn.Flags        = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT;
  ofn.lpstrFile    = spath;
  ofn.lpstrDefExt  = gp.ename;
  ofn.lpstrInitialDir = gp.dname;

  if( GetSaveFileNameA(&ofn) )
  {
    _chdir(init_dir);
    int i=0;
    for( i=0; i<(int)strlen(spath); i++ )
      if( spath[i]== '\\' ) spath[i] = '/';

    strcpy( init_dir, parse_spath(spath).dname );

    for( i=0; i<(int)strlen(init_dir); i++ )
      if( init_dir[i]== '/' ) init_dir[i] = '\\';

    _chdir(curr_dir);

    return true;
  }else
  {
    return false;
  }
}


