using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.Rendering.VirtualTexturing;

public class TemproalVisibility
{
    public Vector3[] ls_v_arr;
    public Vector3[] ls_vn_arr;
    public Vector3[] tn_arr;
    public int[] vtIdx_arr;

    public int[] weld_vtIdx_arr;
    public int[] weld_v_all_arr;
    public int[] weld_v_allInOne_arr;
    public int[] weld_v_perSize_arr;
    public int[] weld_vSize_arr;

    public Vector4[] vfEdges_arr;
    public Vector4[] vfEdgesMap_arr;
    public int[] fEdges_arr;
    public int[] edgeSize_arr;

    public int[] des_ei_arr;
    public float[] des_es_arr;

    public int[] des_vinfo_arr;
    public int[] vis_size_arr;
    public Vector3[] des_vis_single_arr;
    public int[] vis_single_size_arr;

    public bool isInit;
    public bool isGenDone;

    public TemproalVisibility(List<Vector3> ls_v, List<Vector3> ls_vn, List<int> vtIdx, List<Vector3> tn)
    {
        ls_v_arr = ls_v.ToArray();
        ls_vn_arr = ls_vn.ToArray();
        vtIdx_arr = vtIdx.ToArray();
        tn_arr = tn.ToArray();

        weld_vtIdx_arr = new int[vtIdx.Count];
        weld_v_all_arr = new int[vtIdx.Count];
        weld_v_allInOne_arr = new int[vtIdx.Count];
        weld_v_perSize_arr = new int[vtIdx.Count];
        weld_vSize_arr = new int[3];

        vfEdges_arr = new Vector4[1024 * 1024];
        fEdges_arr = new int[vtIdx.Count];
        vfEdgesMap_arr = new Vector4[1024 * 1024];
        edgeSize_arr = new int[1];

        des_ei_arr = new int[1024 * 64];
        des_es_arr = new float[1024 * 64];

        vis_size_arr = new int[2];
        vis_single_size_arr = new int[1];
        des_vis_single_arr = new Vector3[1024 * 1024];

        isInit = false;
        isGenDone = false;
    }

    public void Init()
    {
        GetWeldVertex();
        GetEdgeInfo();

        isInit = true;
    }

    private void GetWeldVertex()
    {
        unsafe
        {
            int nv = ls_v_arr.Length;
            int nf = tn_arr.Length;
            
            fixed (Vector3* ls_vPosPtr = ls_v_arr)                                  // input
            fixed (int* vtIdxPtr = vtIdx_arr)                                       // input
            fixed (int* weld_vtIdxPtr = weld_vtIdx_arr)                             // inout
            fixed (int* weld_vPtr = weld_v_all_arr)                                 // inout
            fixed (int* weld_v_allInOnePtr = weld_v_allInOne_arr)                   // inout
            fixed (int* weld_vPerSizePtr = weld_v_perSize_arr)                      // inout
            fixed (int* weld_vSizePtr = weld_vSize_arr)                             // inout
            {
                weld_func(ls_vPosPtr, vtIdxPtr, nv, nf, weld_vtIdxPtr, weld_vPtr, weld_v_allInOnePtr, weld_vPerSizePtr, weld_vSizePtr);
            }
        }
        Array.Resize(ref weld_v_all_arr, weld_vSize_arr[0]);
        Array.Resize(ref weld_v_allInOne_arr, weld_vSize_arr[1]);
        Array.Resize(ref weld_v_perSize_arr, weld_vSize_arr[2]);
    }

    private void GetEdgeInfo()
    {
        unsafe
        {
            int nv = weld_v_all_arr.Length;
            int nf = tn_arr.Length;

            fixed (int* weld_vtIdxPtr = weld_vtIdx_arr)               // input
            fixed (int* weld_v_allPtr = weld_v_all_arr)               // input
            fixed (Vector4* vfEdgesPtr = vfEdges_arr)                 // output
            fixed (int* fEdgesPtr = fEdges_arr)                       // output
            fixed (Vector4* vfEdgesMapPtr = vfEdgesMap_arr)           // output
            fixed (int* edgeSizePtr = edgeSize_arr)                   // output
            {
                cal_edge_func_ab(weld_vtIdxPtr, weld_v_allPtr, nv, nf, nv, vfEdgesPtr, fEdgesPtr, vfEdgesMapPtr, edgeSizePtr);
            }
        }
        Array.Resize(ref vfEdges_arr, edgeSize_arr[0]);
    }

    public void GetVisibility(int first_vi = 0)
    {
        if (!isInit)
        {
            Debug.LogError("The visibility must be initialized first!");
            return;
        }

        des_vinfo_arr = new int[weld_v_all_arr.Length + 1];
        des_ei_arr = new int[1024 * 64];
        des_es_arr = new float[1024 * 64];

        int nf = tn_arr.Length;
        int dnv = ls_v_arr.Length;
        int nv = weld_v_all_arr.Length;
        int ne = vfEdges_arr.Length;

        unsafe
        {
            fixed (Vector3* ls_vPtr = ls_v_arr)
            fixed (Vector3* ls_vnPtr = ls_vn_arr)
            fixed (int* vtIdxPtr = vtIdx_arr)
            fixed (int* weld_vtIdxPtr = weld_vtIdx_arr)
            fixed (int* weld_v_allPtr = weld_v_all_arr)
            fixed (Vector3* tnPtr = tn_arr)                                    // face normal of non-weld model
            fixed (Vector4* vfEdgesPtr = vfEdges_arr)                          // face normal of non-weld model
            fixed (int* fEdgesPtr = fEdges_arr)

            fixed (int* des_eiPtr = des_ei_arr)                                 // out
            fixed (float* des_esPtr = des_es_arr)                               // out
            fixed (int* des_vinfoPtr = des_vinfo_arr)                           // out
            fixed (int* vis_sizePtr = vis_size_arr)                             // out
            {
                gen_tempVecVis_ondemand_ac(ls_vPtr, ls_vnPtr, vtIdxPtr, weld_vtIdxPtr,
                                           weld_v_allPtr, tnPtr, vfEdgesPtr, vfEdgesPtr, fEdgesPtr,
                                           nf, dnv, nv, ne, first_vi,
                                           des_eiPtr, des_esPtr, des_vinfoPtr, vis_sizePtr);
            }
        }

        Array.Resize(ref des_ei_arr, vis_size_arr[0]);
        Array.Resize(ref des_es_arr, vis_size_arr[0]);

        isGenDone = true;
    }

    public void GetSingleVisibility(int vIdx)
    {
        if (!isGenDone)
        {
            Debug.LogError("Please generate the completed visibility by calling GenVisibility() " +
                "           before inquiring the visibility of a single vertex.");
            return;
        }
        if (!isInit)
        {
            Debug.LogError("The visibility must be initialized first!");
            return;
        }

        int nv = weld_v_all_arr.Length;    
        int nf = tn_arr.Length;
        int nei = des_ei_arr.Length;
        int dnv = ls_v_arr.Length;            
        int ne = vfEdgesMap_arr.Length;

        des_vis_single_arr = new Vector3[1024 * 1024];

        unsafe
        {
            fixed (Vector3* ls_vPtr = ls_v_arr)
            fixed (Vector3* ls_vnPtr = ls_vn_arr)
            fixed (int* weld_v_allPtr = weld_v_all_arr)
            fixed (Vector4* vfEdgesPtr = vfEdges_arr)
            fixed (Vector4* vfEdgesMapPtr = vfEdgesMap_arr)
            fixed (int* des_eiPtr = des_ei_arr)
            fixed (float* des_esPtr = des_es_arr)
            fixed (int* des_vinfoPtr = des_vinfo_arr)
            fixed (Vector3* des_vis_singlePtr = des_vis_single_arr)               // output
            fixed (int* vis_single_sizePtr = vis_single_size_arr)                 // output
            {
                get_singleV_vecvis_ac(ls_vPtr, ls_vnPtr, weld_v_allPtr, vfEdgesMapPtr, des_eiPtr, des_esPtr, des_vinfoPtr,
                                      nv, nf, nei, dnv, ne,
                                      vIdx,
                                      des_vis_singlePtr, vis_single_sizePtr);
            }
        }
        Array.Resize(ref des_vis_single_arr, vis_single_size_arr[0]);
    }

    #region vulkan version native plugin
    [DllImport("Vulkan_temvecvis_gen", EntryPoint = "weld_func")]
    private unsafe static extern void weld_func(Vector3* _ls_v,                      // vpos 
                                                int* _vtIdx,
                                                int _dnv, int _nf,                   // _dnv: non-welded vertex number
                                                int* _weld_vtIdx,                    // output, length = 1728
                                                int* _weld_vidx_all,                 // output, length = 1728
                                                int* _weld_vidx_allInOne,            // output, length = 1728
                                                int* _weld_vidx_perSize,             // output, length = 1728
                                                int* _weld_vSize);                   // output, length = 3

    [DllImport("Vulkan_temvecvis_gen", EntryPoint = "cal_edge_func_ab")]
    private unsafe static extern void cal_edge_func_ab(int* weld_vtidx,                      // vtidx after welding
                                                        int* _weld_vidx_all,
                                                        int _nv, int _nf, int _weld_nv,      // _nv: welded vertex number = 290, _nf = #triangle
                                                        Vector4* _vfEdges,                   // output, fixed size = 1024*1024 = 1048576
                                                        int* _fEdges,                        // output, length = 1728
                                                        Vector4* _vfEdgesMap,                // output, fixed size = 1024*1024 = 1048576
                                                        int* _edgeSize);                     // output, length = 1

    [DllImport("Vulkan_temvecvis_gen", EntryPoint = "gen_tempVecVis_ondemand_ac")]
    private unsafe static extern void gen_tempVecVis_ondemand_ac(Vector3* _ls_v,
                                                                 Vector3* _ls_vn,
                                                                 int* _vtIdx,
                                                                 int* _weld_vtIdx,
                                                                 int* _weld_vidx_all,
                                                                 Vector3* _tNormal,
                                                                 Vector4* _vfEdges,
                                                                 Vector4* _vfEdgesMap,
                                                                 int* _fEdges,
                                                                 int _nf, int _dnv,          // _nf = #triangle, _dnv: non-welded vertex number
                                                                 int _nv, int _ne,           // _nv: welded vertex number = 290, _ne: _vfEdges.size
                                                                 int _first_vi,           // _nv: welded vertex number = 290, _ne: _vfEdges.size
                                                                 int* des_ei,                // output, fixed size = 1024*1024*nv = 304087040
                                                                 float* des_es,                // output, fixed size = 1024*1024*nv = 304087040
                                                                 int* des_vinfo0,            // output, length = nv+1
                                                                 int* vis_size);             // output, length = 1    


    [DllImport("Vulkan_temvecvis_gen", EntryPoint = "get_singleV_vecvis_ac")]
    private unsafe static extern void get_singleV_vecvis_ac(Vector3* _ls_v,             // vpos 
                                                            Vector3* _ls_vn,            // vnormal
                                                            int* _weld_vidx_all,        // 290
                                                            Vector4* _vfEdgesMap,       // 864
                                                            int* des_ei,                // 172838 
                                                            float* des_es,              // 172838
                                                            int* des_vinfo0,            // length = nv+1 = 291
                                                            int _nv, int _nf,           // _nv: welded vertex number = 290, _nf = #triangle = 576
                                                            int _nei, int _dnv,         // _nei = 172838,  _dnv: non-welded vertex number = 1728
                                                            int _neMap,                 // _ne: _vfEdges.size
                                                            int _vIdx,                  // current vertex idx of visibility
                                                            Vector3* des_vis_single,    // output: fixed length = 1024*1024
                                                            int* vis_single_size);      // output: fixed length = 1
    #endregion
}
