using UnityEngine;
using System.IO;
using System.Collections.Generic;
using UnityEngine.UIElements;
using System.Linq;

namespace q_common
{
    public static class Model
    {
        public static void Unitize(ref List<Vector3>vertex, Vector3 center, Vector3 max, Vector3 min)
        {
            Vector3 distance = max - min;
            float scale = 2 / Mathf.Max( distance.x, Mathf.Max( distance.y, distance.z ) );

            for (int i = 0; i < vertex.Count; i++)
            { 
                vertex[i] -= center;
                vertex[i] *= scale;
            }
        }

        public static Vector3 Unitize(Vector3 vertex, Vector3 center, Vector3 max, Vector3 min)
        {
            Vector3 distance = max - min;
            float scale = 2 / Mathf.Max(distance.x, Mathf.Max(distance.y, distance.z));

            vertex -= center;
            vertex *= scale;

            return vertex;
        }

        public static Transform Unitize(Transform joints, Vector3 center, Vector3 max, Vector3 min)
        {
            Debug.Log(joints.transform.position);
            Debug.Log(joints.transform.localScale);
            Vector3 distance = max - min;
            float scale = 2 / Mathf.Max(distance.x, Mathf.Max(distance.y, distance.z));

            //joints.position -= center;
            joints.localScale *= scale;

            Debug.Log(joints.transform.position);
            Debug.Log(joints.transform.localScale);

            return joints;
        }

        public static void GetBounds(List<Vector3> vertex, ref Vector3 max, ref Vector3 min, ref Vector3 center)
        {
            max = vertex[0];
            min = vertex[0];

            // Find the max and min
            for (int i = 1; i < vertex.Count; i++)
            {
                max.x = Mathf.Max( vertex[i].x, max.x );
                max.y = Mathf.Max( vertex[i].y, max.y );
                max.z = Mathf.Max( vertex[i].z, max.z );

                min.x = Mathf.Min(vertex[i].x, min.x);
                min.y = Mathf.Min(vertex[i].y, min.y);
                min.z = Mathf.Min(vertex[i].z, min.z);
            }

            //Debug.Log(max + "," + min);
            center = (max + min) / 2.0f;

            Vector3 distance = max - min;
            float scale = 2 / Mathf.Max(distance.x, Mathf.Max(distance.y, distance.z));
        }

        public static void GetBoundsWS(List<Vector3> vertex, ref Vector3 max, ref Vector3 min, ref Vector3 center)
        {
            max = vertex[0];
            min = vertex[0];

            // Find the max and min
            for (int i = 1; i < vertex.Count; i++)
            {
                max.x = Mathf.Abs(vertex[i].x) >= Mathf.Abs(max.x) ? vertex[i].x : max.x;
                max.y = Mathf.Abs(vertex[i].y) >= Mathf.Abs(max.y) ? vertex[i].y : max.y;
                max.z = Mathf.Abs(vertex[i].z) >= Mathf.Abs(max.z) ? vertex[i].z : max.z;

                min.x = Mathf.Abs(vertex[i].x) <= Mathf.Abs(min.x) ? vertex[i].x : min.x;
                min.y = Mathf.Abs(vertex[i].y) <= Mathf.Abs(min.y) ? vertex[i].y : min.y;
                min.z = Mathf.Abs(vertex[i].z) <= Mathf.Abs(min.z) ? vertex[i].z : min.z;
            }

            Debug.Log(max + "," + min);
            center = (max + min) / 2.0f;

            Vector3 distance = max - min;
            float scale = 2 / Mathf.Max(distance.x, Mathf.Max(distance.y, distance.z));
        }
    }
}
