using System;
using UnityEngine;

[Serializable]
public class qVtAnimData : ScriptableObject
{
    [Serializable]
    public class qFrameData
    {
        public float time;
        public Vector3[] animVertex;
        public Vector3[] animVNormal;
    }

    [SerializeField]
    public string animName;
    [SerializeField]
    public float animLen;
    [SerializeField]
    public int frame;
    [SerializeField]
    public qFrameData[] frameDatas;
}