using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.UI;
using UnityEditor;
using static qVtAnimData;

#if UNITY_EDITOR
public class ConvertAnimEditor : EditorWindow
{
    private qAnimDataAssets _AnimDataAssets;
    private HashSet<string> _selectSet;
    private Mesh _mesh;

    [MenuItem("Tools/ConvertAnimEditor")]
    public static void Open()
    {
        GetWindow<ConvertAnimEditor>().Show();
    }

    private void OnGUI()
    {
        _AnimDataAssets = (qAnimDataAssets)EditorGUILayout.ObjectField("qObjectAnimationAsset: ", _AnimDataAssets, typeof(qAnimDataAssets), false);

        if (_AnimDataAssets == null)
        {
            EditorGUILayout.HelpBox("Please select a prefab containing the component of Mesh Filter, Mesh Renderer and AnimDataAssets", MessageType.Error);
            return;
        }
        if (_AnimDataAssets.gameObject.GetComponent<MeshFilter>() == null)
        {
            EditorGUILayout.HelpBox("Please select a prefab containing the component of Mesh Filter!", MessageType.Error);
            return;
        }
        else
        {
            //_mesh = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh;
            _mesh = new();
            _mesh.vertices = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh.vertices;
            _mesh.normals = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh.normals;
            _mesh.bindposes = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh.bindposes;
            _mesh.boneWeights = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh.boneWeights;
            _mesh.triangles = _AnimDataAssets.gameObject.GetComponent<MeshFilter>().sharedMesh.triangles;
        }

        GUILayout.Space(5);
        GUILayout.Label("AnimDataAssets List:");
        int index = 0;

        foreach (AnimData a in _AnimDataAssets.animDatasList)
        {
            EditorGUILayout.BeginHorizontal();
            {
                bool isSelect = _selectSet.Contains(a.animName);
                bool newState = EditorGUILayout.Toggle(isSelect);
                if (isSelect != newState)
                {
                    if (!newState)
                    {
                        _selectSet.Remove(a.animName);
                    }
                    else
                    {
                        _selectSet.Add(a.animName);
                    }
                }
                EditorGUILayout.ObjectField(a, typeof(AnimData), false);
            }
            EditorGUILayout.EndHorizontal();
            index++;
        }
        GUILayout.Space(5);
        if (index == 0)
        {
            EditorGUILayout.HelpBox("AnimDataAssets List is null!", MessageType.Error);
            return;
        }
        GUILayout.Space(5);
        if (GUILayout.Button("Convert to vertex animation"))
        {
            string dir = EditorUtility.SaveFolderPanel("Convert to vertex animation", "", "");
            if (!string.IsNullOrEmpty(dir))
            {
                dir = dir.Replace("\\", "/");
                if (!dir.StartsWith(Application.dataPath))
                {
                    return;
                }
                dir = dir.Replace(Application.dataPath, "Assets");
                ExportAnim(dir);
            }
        }
    }

    private void OnEnable()
    {
        _mesh = null;
        _selectSet = new HashSet<string>();
    }

    private void ExportAnim(string dir)
    {
        foreach (AnimData animData in _AnimDataAssets.animDatasList)
        {
            if (_selectSet.Contains(animData.animName))
                ExportAnimData(animData, dir);
        }
    }

    private void ExportAnimData(AnimData m_animData, string dir)
    {
        qVtAnimData vtAnimData = ScriptableObject.CreateInstance<qVtAnimData>();
        vtAnimData.animName = m_animData.animName;
        vtAnimData.animLen = m_animData.animLen;
        vtAnimData.frame = m_animData.frame;
        vtAnimData.frameDatas = Skinned2VertexAnim(m_animData);

        string path = dir + "/" + m_animData.animName + ".asset";
        AssetDatabase.CreateAsset(vtAnimData, path);
    }

    public qVtAnimData.qFrameData[] Skinned2VertexAnim(AnimData m_animData)
    {
        
        List<qVtAnimData.qFrameData> qFrameDataList = new();
        for (int i = 0; i < m_animData.frameDatas.Length; i++)
        {
            qVtAnimData.qFrameData m_qFrameData = new();
            m_qFrameData.time = m_animData.frameDatas[i].time;
            m_qFrameData.animVertex = GetVertexPosFromAnims(i, m_animData, _mesh);

            Mesh _tempMesh = new();
            _tempMesh.vertices = m_qFrameData.animVertex;
            _tempMesh.normals = _mesh.normals;
            _tempMesh.triangles = _mesh.triangles;
            _tempMesh.RecalculateNormals();

            m_qFrameData.animVNormal = _tempMesh.normals;

            qFrameDataList.Add(m_qFrameData);
        }

        return qFrameDataList.ToArray();
    }

    Vector3[] GetVertexPosFromAnims(int _frame, AnimData animData, Mesh _mesh)
    {
        Vector3[] _srcPoints = _mesh.vertices;
        Matrix4x4[] _bindPoses = _mesh.bindposes;
        Vector3[] _tempPoints = new Vector3[_srcPoints.Length];

        AnimData.FrameData frameData = animData.frameDatas[_frame];
        for (int i = 0; i < _srcPoints.Length; ++i)
        {
            Vector3 point = _srcPoints[i];
            BoneWeight weight = _mesh.boneWeights[i];
            Matrix4x4 tempMat0 = frameData.matrix4X4s[weight.boneIndex0] * _bindPoses[weight.boneIndex0];
            Matrix4x4 tempMat1 = frameData.matrix4X4s[weight.boneIndex1] * _bindPoses[weight.boneIndex1];
            Matrix4x4 tempMat2 = frameData.matrix4X4s[weight.boneIndex2] * _bindPoses[weight.boneIndex2];
            Matrix4x4 tempMat3 = frameData.matrix4X4s[weight.boneIndex3] * _bindPoses[weight.boneIndex3];

            Vector3 temp = tempMat0.MultiplyPoint(point) * weight.weight0 +
                                   tempMat1.MultiplyPoint(point) * weight.weight1 +
                                   tempMat2.MultiplyPoint(point) * weight.weight2 +
                                   tempMat3.MultiplyPoint(point) * weight.weight3;

            _tempPoints[i] = temp;
        }
        Debug.Log(_tempPoints.Length);

        return _tempPoints;
    }
}
#endif