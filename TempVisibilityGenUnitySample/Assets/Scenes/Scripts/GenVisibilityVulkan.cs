using System;
using System.Collections.Generic;
using System.Linq;
using TMPro;
using UnityEngine;

public class GenVisibilityVulkan : MonoBehaviour
{
    #region public variable
    [Header("Model")]
    public GameObject[] qObjects;

    [Header("Visibility Mode")]
    [NonSerialized] public bool isDrawVisibility;
    public GameObject lineRendererContainer;
    public GameObject vertexRenderer;
    public GameObject lineRendererPrefab;
    [Range(0, 407)]
    [NonSerialized] public int shownVIdx = 331;

    [Header("Depth View")]
    [Tooltip("If the camera will transform to the shading point view.")]
    public bool isDepthCamera = false;
    public Camera depthCamera;
    #endregion

    #region variable declaration
    List<Vector3> ls_v = new List<Vector3>();
    List<Vector3> ls_vn = new List<Vector3>();
    List<int> vtIdx = new List<int>();
    List<Vector3> tn = new List<Vector3>();

    qObjectManager m_qObjectManager;
    [NonSerialized] public TemproalVisibility tv;

    // Animation
    [NonSerialized] public bool isPlayingAnimation;
    [NonSerialized] public bool playNextFrame;
    [NonSerialized] public int AnimIdx = 2;
    [NonSerialized] public int currentFrame;
    [NonSerialized] public List<qVtAnimData.qFrameData[]> frameDatasWS;

    int lastShownVIdx;
    int lastFrame;
    List<qVtAnimData> objAnimOSList;
    List<List<Vector3[]>> frameTNormalWS;
    GameObject q_obj;
    MeshFilter q_filter;
    List<Vector3> animTempVPosOS;
    float _time;
    int frameCount;
    #endregion

    void Start()
    {
        // Initialize the vertex animation
        frameDatasWS = new();
        frameTNormalWS = new();
        objAnimOSList = new List<qVtAnimData>();
        LoadObjVertexAnimation(ref objAnimOSList);

        // Set up camera and default state
        isDrawVisibility = true;
        frameCount = 0;
        currentFrame = -1;
        _time = 0;
        lastFrame = currentFrame;
        lastShownVIdx = shownVIdx;
        isPlayingAnimation = false;
        playNextFrame = false;

        // merge all objects into one
        GameObject q_obj = new GameObject("Merge Object");
        q_obj.AddComponent<MeshFilter>();
        q_obj.AddComponent<MeshRenderer>();
        q_filter = q_obj.GetComponent<MeshFilter>();

        // unitize object 
        m_qObjectManager = new qObjectManager();
        m_qObjectManager.UnitizeObjects(ref frameDatasWS, qObjects, q_obj, objAnimOSList);

        // hide the original obejct and reamin the merge one
        foreach (var obj in qObjects)
            obj.SetActive(false);
        
        // prepare the data to calculate visibility
        ls_v = m_qObjectManager.vPos_unitize;
        ls_vn = m_qObjectManager.vNormals;
        vtIdx = m_qObjectManager.vtIdx;
        tn = m_qObjectManager.GetFaceNormal(ls_v.ToArray(), vtIdx.ToArray(), m_qObjectManager.nTriangle).ToList();

        for (int i = 0; i < frameDatasWS.Count; i++)
        {
            int frameLen = frameDatasWS[i].Length;
            List<Vector3[]> m_tNormalList = new();
            for (int k = 0; k < frameLen; k++)
            {
                Vector3[] m_tn = m_qObjectManager.GetFaceNormal(frameDatasWS[i][k].animVertex, vtIdx.ToArray(), m_qObjectManager.nTriangle);
                m_tNormalList.Add(m_tn);
            }
            frameTNormalWS.Add(m_tNormalList);
        }

        #region visibility
        tv= new(ls_v, ls_vn, vtIdx, tn);

        tv.Init();
        tv.GetVisibility();
        if (isDrawVisibility)
            tv.GetSingleVisibility(shownVIdx);
        #endregion

        if (isDrawVisibility)
            DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
        UpdateCamera(shownVIdx);
    }

    void Update()
    {
        if (lastShownVIdx != shownVIdx)
        {
            if (isDrawVisibility) 
            { 
                tv.GetSingleVisibility(shownVIdx);
                DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
            }
            UpdateCamera(shownVIdx);
            lastShownVIdx = shownVIdx;
        }

        // Debug manually
        if (playNextFrame)
        {
            GameObject m_object = GameObject.Find("Merge Object");
            if (m_object.GetComponent<MeshCollider>() != null)
                Destroy(m_object.GetComponent<MeshCollider>());

            if (currentFrame < frameDatasWS[AnimIdx].Length - 1)
                currentFrame++;
            else
                currentFrame = 0;
            tv.ls_v_arr = frameDatasWS[AnimIdx][currentFrame].animVertex;
            tv.ls_vn_arr = frameDatasWS[AnimIdx][currentFrame].animVNormal;
            tv.tn_arr = frameTNormalWS[AnimIdx][currentFrame];
            Play(currentFrame);
            m_object.AddComponent<MeshCollider>();
            if (frameCount % 8 == 0)
                tv.GetVisibility();
            else if (frameCount % 2 == 0)
                tv.GetVisibility(tv.weld_v_all_arr.Length - 121);

            if (isDrawVisibility)
            {
                tv.GetSingleVisibility(shownVIdx);
                DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
            }

            //Debug.Log("Frame:" + currentFrame);
            lastFrame = currentFrame;
            frameCount++;
            playNextFrame = false;
        }


        if (isPlayingAnimation)
        {
            _time += Time.deltaTime;
            _time %= objAnimOSList[AnimIdx].animLen;
            currentFrame = (int)(_time / (1.0f / objAnimOSList[AnimIdx].frame));
            Play(currentFrame);
            tv.ls_v_arr = frameDatasWS[AnimIdx][currentFrame].animVertex;
            tv.ls_vn_arr = frameDatasWS[AnimIdx][currentFrame].animVNormal;
            tv.tn_arr = frameTNormalWS[AnimIdx][currentFrame];

            if (frameCount % 8 == 0)
                tv.GetVisibility();
            else if (frameCount % 2 == 0)
                tv.GetVisibility(shownVIdx);


            if (isDrawVisibility)
            {
                tv.GetSingleVisibility(shownVIdx);
                DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
            }
            UpdateCamera(shownVIdx);

            frameCount++;
            lastFrame = currentFrame;
        }

        // restore the default state
        if (Input.GetKeyDown(KeyCode.F3))
            ResetModel();
    }

    void Play(int frame)
    {
        q_filter.sharedMesh.SetVertices(frameDatasWS[AnimIdx][frame].animVertex);
        q_filter.sharedMesh.SetNormals(frameDatasWS[AnimIdx][frame].animVNormal);
    }

    void DrawVisibilityPairs(int m_shownVIdx, List<Vector3> visibilityArr)
    {
        GameObject lineObj;
        LineRenderer m_lineRenderer;

        // Remove all the child object of lineRendererContainer
        if (lastShownVIdx != shownVIdx)
        {
            Transform m_transform;
            for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
            {
                m_transform = lineRendererContainer.transform.GetChild(i);
                GameObject.Destroy(m_transform.gameObject);
            }
        }

        if (lastFrame != currentFrame)
        {
            Transform m_transform;
            for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
            {
                m_transform = lineRendererContainer.transform.GetChild(i);
                GameObject.Destroy(m_transform.gameObject);
            }
        }

        for (int i = 0; i < visibilityArr.Count; i+=2)
        {
            lineObj = Instantiate(lineRendererPrefab);
            lineObj.transform.parent = lineRendererContainer.transform;
            m_lineRenderer = lineObj.GetComponent<LineRenderer>();
            m_lineRenderer.positionCount = 2;
            m_lineRenderer.SetPosition(0, visibilityArr[i]);
            m_lineRenderer.SetPosition(1, visibilityArr[i+1]);
        } 

        // draw the position
        vertexRenderer.transform.position = tv.ls_v_arr[tv.weld_v_all_arr[m_shownVIdx]];
    }

    void UpdateCamera(int m_shownVIdx)
    {
        depthCamera.transform.position = tv.ls_v_arr[tv.weld_v_all_arr[m_shownVIdx]];
        depthCamera.transform.LookAt(new Vector3(0, 0, 0));
    }

    public void ResetModel()
    {
        isPlayingAnimation = false;

        Transform m_transform;
        for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
        {
            m_transform = lineRendererContainer.transform.GetChild(i);
            GameObject.Destroy(m_transform.gameObject);
        }

        q_filter.sharedMesh.SetVertices(ls_v);

        tv.ls_v_arr = ls_v.ToArray();
        tv.ls_vn_arr = ls_vn.ToArray();
        tv.tn_arr = tn.ToArray();
        tv.vtIdx_arr = vtIdx.ToArray();

        tv.GetVisibility();
        if (isDrawVisibility)
        { 
            tv.GetSingleVisibility(shownVIdx);
            DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
        }
        UpdateCamera(shownVIdx);
    }

    public void OnShowFrame()
    {
        if (currentFrame != -1)
        {
            Play(currentFrame);

            tv.ls_v_arr = frameDatasWS[AnimIdx][currentFrame].animVertex;
            tv.ls_vn_arr = frameDatasWS[AnimIdx][currentFrame].animVNormal;
            tv.tn_arr = frameTNormalWS[AnimIdx][currentFrame];

            if (frameCount % 8 == 0)
                tv.GetVisibility();
            else if (frameCount % 2 == 0)
                tv.GetVisibility(shownVIdx);


            if (isDrawVisibility)
            {
                tv.GetSingleVisibility(shownVIdx);
                DrawVisibilityPairs(shownVIdx, tv.des_vis_single_arr.ToList());
            }
            UpdateCamera(shownVIdx);
        }
        else
            ResetModel();
    }

    public void LoadObjVertexAnimation(ref List<qVtAnimData> objVtAnimDataList)
    {
        foreach (GameObject qObj in qObjects)
        {
            if (qObj.gameObject.GetComponent<qObjectAnimationAsset>() != null)
                objVtAnimDataList = qObj.gameObject.GetComponent<qObjectAnimationAsset>().vtAnimDatasList;

        }
    }

    public void UpdateVisibilityWidth(float width)
    {
        Transform m_line;
        LineRenderer m_lineRenderer;
        for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
        {
            m_line = lineRendererContainer.transform.GetChild(i);
            m_lineRenderer = m_line.GetComponent<LineRenderer>();
            m_lineRenderer.startWidth = width;
            m_lineRenderer.endWidth = width;
        }
    }

    public void UpdateVisibilityColor(Color color)
    {
        Transform m_line;
        LineRenderer m_lineRenderer;
        for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
        {
            m_line = lineRendererContainer.transform.GetChild(i);
            m_lineRenderer = m_line.GetComponent<LineRenderer>();
            m_lineRenderer.startColor = color;
            m_lineRenderer.endColor = color;
        }
    }

    public void DestroyVisibilityDraw()
    {
        Transform m_transform;
        for (int i = 0; i < lineRendererContainer.transform.childCount; i++)
        {
            m_transform = lineRendererContainer.transform.GetChild(i);
            GameObject.Destroy(m_transform.gameObject);
        }
    }


}