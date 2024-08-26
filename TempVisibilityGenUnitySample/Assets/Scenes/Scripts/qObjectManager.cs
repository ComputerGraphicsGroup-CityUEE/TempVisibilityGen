using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using Unity.VisualScripting;
using UnityEngine.Analytics;


public class qObjectManager
{
    public List<Vector3> vPos;
    public List<Vector3> vNormals;
    public List<Vector3> tNormals;
    public List<Vector4> edges;
    public List<Vector3> vPos_world;
    public List<Vector3> vPos_unitize;
    public List<int> vtIdx;
    public List<int> tIdx;
    public List<int> vIdx;

    public Vector3[] vPos_arr;
    public Vector3[] vNormals_arr;
    public Vector3[] tNormals_arr;
    public Vector4[] edges_arr;
    public Vector4[] vPos_world_arr;
    public Vector4[] vPos_unitize_arr;
    public int[] vtIdx_arr;
    public int[] tIdx_arr;
    public int[] vIdx_arr;

    public int nTriangle;

    public Vector3 boundMax;
    public Vector3 boundMin;
    public Vector3 objCenter;

    List<int> vPosLastIdxByObj = new List<int>();
    List<int> triangleLastIdxByObj = new List<int>();

    public Matrix4x4 localToWorld;
    public Matrix4x4 worldToLocal;

    Matrix4x4[] bindPoses;
    Mesh mesh;

    public qObjectManager()
    {
        vPos = new List<Vector3>();
        vNormals = new List<Vector3>();
        tNormals = new List<Vector3>();
        edges = new List<Vector4>();
        vPos_world = new List<Vector3>();
        vPos_unitize = new List<Vector3>();
        vtIdx = new List<int>();
        tIdx = new List<int>();
        vIdx = new List<int>();

        nTriangle = 0;
    }

    public void UnitizeObjects(ref List<qVtAnimData.qFrameData[]> frameDatasWS_unitized,
                                GameObject[] objects, GameObject merge_obj, List<qVtAnimData> objAnimOSList)
    {
        vPos.Clear();
        vNormals.Clear();
        vtIdx.Clear();
        tIdx.Clear();
        int objIdx = 0;

        List<Vector3> vPos_plane = new List<Vector3>();
        List<Vector3> vPos_wolrd_plane = new List<Vector3>();
        List<Vector3> vNormals_plane = new List<Vector3>();

        int meshFiltersN = objects.Length;
        CombineInstance[] combine = new CombineInstance[meshFiltersN];
        List<Material> materials = new List<Material>();

        /////////////// Prepare Static Model Data ///////////////
        foreach (GameObject obj in objects)
        {
            if (obj.GetComponent<MeshFilter>() != null)
            {
                //Debug.Log("MeshRenderer was found on " + obj.gameObject.name);
                mesh = obj.GetComponent<MeshFilter>().sharedMesh;
            }
            else
            {
                //Debug.Log("No MeshRenderer was found on " + obj.gameObject.name);
                return;
            }

            localToWorld = obj.transform.localToWorldMatrix;
            bindPoses = obj.GetComponent<MeshFilter>().mesh.bindposes;
            boundMax = new Vector3(mesh.bounds.max.x, mesh.bounds.max.y, mesh.bounds.max.z);
            boundMin = new Vector3(mesh.bounds.min.x, mesh.bounds.min.y, mesh.bounds.min.z);

            int firstVertex = vPos.Count;
            int firstTriangle = tIdx.Count;

            //Debug.Log(firstVertex);
            var tempT = Enumerable.Range(0, mesh.triangles.Count() / 3).ToArray();      // tempT.count = how many triangle we have got
            nTriangle += tempT.Length;

            vPos.AddRange(mesh.vertices);
            vPos_world.AddRange(mesh.vertices.Select(vec => localToWorld.MultiplyPoint3x4(vec)));

            if (!IsAnimModel(obj))
            {
                vPos_plane.AddRange(mesh.vertices);
                vPos_wolrd_plane.AddRange(mesh.vertices.Select(vec => localToWorld.MultiplyPoint3x4(vec)));
                vNormals_plane.AddRange(mesh.normals);
            }

            vNormals.AddRange(mesh.normals);

            tIdx.AddRange(tempT.Select(index => index + firstTriangle));
            vtIdx.AddRange(mesh.triangles.Select(index => index + firstVertex));

            vPosLastIdxByObj.Add(vPos.Count);
            triangleLastIdxByObj.Add(tIdx.Count);
            objIdx++;
        }

        /////////////// Unitize the model ////////////////////
        q_common.Model.GetBounds(vPos_world, ref boundMax, ref boundMin, ref objCenter);
        //vPos_unitize.AddRange(vPos_world);
        vPos_unitize.AddRange(vPos_world.Select(vec => q_common.Model.Unitize(vec, objCenter, boundMax, boundMin)));

        /////////////// Save the unitized vertex os back to the mesh //////
        int i = 0;
        foreach (GameObject obj in objects)
        {
            if (obj.GetComponent<MeshFilter>() != null)
                mesh = obj.GetComponent<MeshFilter>().sharedMesh;
            else
                mesh = obj.GetComponent<SkinnedMeshRenderer>().sharedMesh;

            MeshRenderer renders = obj.GetComponent<MeshRenderer>();

            List<Vector3> m_vPos_unitizedOS = new List<Vector3>();
            List<Vector3> m_vPos_unitizedWS = new List<Vector3>();

            if (i == 0)
            {
                m_vPos_unitizedWS = vPos_unitize.GetRange(0, vPosLastIdxByObj[i]);
            }
            else
            {
                m_vPos_unitizedWS = vPos_unitize.GetRange(vPosLastIdxByObj[i - 1], vPosLastIdxByObj[i] - vPosLastIdxByObj[i - 1]);
            }

            worldToLocal = obj.transform.worldToLocalMatrix;
            m_vPos_unitizedOS.AddRange(m_vPos_unitizedWS.Select(vec => worldToLocal.MultiplyPoint3x4(vec)));

            mesh.vertices = m_vPos_unitizedOS.ToArray();
            mesh.RecalculateBounds();

            ///////////////////////// Combine them into a single mesh /////////////////////////////////

            combine[i].mesh = mesh;
            combine[i].transform = obj.transform.localToWorldMatrix;
            materials.Add(renders.sharedMaterial);

            i++;
        }

        /////////////// Combine obj and plane as a single mesh ///////////
        {
            Mesh combine_mesh = new Mesh();
            combine_mesh.CombineMeshes(combine, false, true);

            MeshFilter q_filter = merge_obj.GetComponent<MeshFilter>();
            MeshRenderer q_meshrenderer = merge_obj.GetComponent<MeshRenderer>();
            q_filter.sharedMesh = combine_mesh;
            q_meshrenderer.materials = materials.ToArray();
        }

        /////////////// Prepare Animation ////////////////////
        List<qVtAnimData.qFrameData[]> frameDatasWS = new();
        int animCount = objAnimOSList.Count;
        foreach (GameObject obj in objects)
        {
            if (obj.GetComponent<MeshFilter>() != null)
            {
                //Debug.Log("MeshRenderer was found on " + obj.gameObject.name);
                mesh = obj.GetComponent<MeshFilter>().sharedMesh;
            }
            else
            {
                //Debug.Log("No MeshRenderer was found on " + obj.gameObject.name);
                return;
            }

            localToWorld = obj.transform.localToWorldMatrix;

            if (IsAnimModel(obj))
            {
                for (int j = 0; j < animCount; j++)
                {
                    int frameLen = objAnimOSList[j].frameDatas.Length;
                    List<qVtAnimData.qFrameData> m_qFrameDataList = new();
                    for (int k = 0; k < frameLen; k++)
                    {
                        List<Vector3> temp_animVertex = new();
                        List<Vector3> temp_animVNormal = new();
                        qVtAnimData.qFrameData src_frameData = objAnimOSList[j].frameDatas[k];
                        qVtAnimData.qFrameData des_frameData = new();

                        temp_animVertex.AddRange(src_frameData.animVertex.Select(vec => localToWorld.MultiplyPoint3x4(vec)));
                        temp_animVertex.AddRange(vPos_wolrd_plane);

                        temp_animVNormal.AddRange(src_frameData.animVNormal);
                        temp_animVNormal.AddRange(vNormals_plane);

                        des_frameData.time = src_frameData.time;
                        des_frameData.animVertex = temp_animVertex.ToArray();
                        des_frameData.animVNormal = temp_animVNormal.ToArray();

                        m_qFrameDataList.Add(des_frameData);
                    }
                    frameDatasWS.Add(m_qFrameDataList.ToArray());
                }
            }
        }

        // Debug.Log(boundMax + ", " + boundMin + ", " + objCenter);

        /////////////// Unitize the Animation /////////////////
        for (int j = 0; j < animCount; j++)
        {
            int frameLen = frameDatasWS[j].Length;
            List<qVtAnimData.qFrameData> m_qFrameDataList = new();
            for (int k = 0; k < frameLen; k++)
            {
                List<Vector3> temp_animVertex = new();
                List<Vector3> temp_animVNormal = new();
                qVtAnimData.qFrameData src_frameData = frameDatasWS[j][k];
                qVtAnimData.qFrameData des_frameData = new();

                temp_animVertex.AddRange(src_frameData.animVertex.Select(vec => q_common.Model.Unitize(vec, objCenter, boundMax, boundMin)));
                //temp_animVertex.AddRange(src_frameData.animVertex);
                temp_animVNormal.AddRange(src_frameData.animVNormal);
                des_frameData.time = src_frameData.time;
                des_frameData.animVertex = temp_animVertex.ToArray();
                des_frameData.animVNormal = temp_animVNormal.ToArray();

                m_qFrameDataList.Add(des_frameData);
            }
            frameDatasWS_unitized.Add(m_qFrameDataList.ToArray());
        }
    }

    public bool IsAnimModel(GameObject obj)
    {
        if (obj.GetComponent<qObjectAnimationAsset>() != null)
            return true;
        else
            return false;
    }


    public Vector3[] GetFaceNormal(Vector3[] m_vPos, int[] m_vtIdx, int m_nTriangle)
    {
        List<Vector3> m_tn = new();
        for (int i = 0; i < m_nTriangle; i++)
        {
            m_tn.Add(Vector3.Cross(m_vPos[m_vtIdx[i * 3 + 1]] - m_vPos[m_vtIdx[i * 3 + 0]],
                                   m_vPos[m_vtIdx[i * 3 + 2]] - m_vPos[m_vtIdx[i * 3 + 1]]).normalized);
        }

        return m_tn.ToArray();
    }
}