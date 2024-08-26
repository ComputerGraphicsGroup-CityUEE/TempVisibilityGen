#ifndef G_FILEDIALOG_H
#define G_FILEDIALOG_H

#define G_FILEDIALOG_NMAXFILE 256

class GFileDialog
{
  public:

    char spath[256];
    char init_dir[256];
    char spath_relative[256];

    GFileDialog();

    bool prompt_open( const char *file_filter );
    bool prompt_save( const char *file_filter );

};



//char G_FILEDIALOG_ALLFILE[] = "All files\0*.*\0\0";
//bool g_dialog_open( char *spath, const char *file_filter=G_FILEDIALOG_ALLFILE );


#endif
