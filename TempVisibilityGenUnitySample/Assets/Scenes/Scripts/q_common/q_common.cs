using System;
using System.IO;
using TMPro;
using UnityEngine;

namespace q_common
{
    public static class q_common
    {
        public static string ReadLine(BinaryReader reader)
        {
            char nextChar;
            string line = "";

            while ((nextChar = reader.ReadChar()) != '\n')
            {
                line += nextChar;
            }

            return line;
        }
        public static float ReverseBytes(float value)
        {
            byte[] bytes = BitConverter.GetBytes(value);
            Array.Reverse(bytes);
            return BitConverter.ToSingle(bytes, 0);
        }

        public static Vector3 mix(Vector3 a, Vector3 b, float z)
        {
            return new Vector3(a.x * (1 - z) + b.x * z, a.y * (1 - z) + b.y * z, a.z * (1 - z) + b.z * z);
        }

        public static void get_coordinate_system(Vector3 ay, ref Vector3 bx, ref Vector3 by, ref Vector3 bz)
        {
            by = Vector3.Normalize(ay);
            bx = new Vector3(1, 0, 0);
            Vector3 temp = new Vector3(by.z, 0.0f, -by.x);
            Vector3 refere = mix(temp, bx, Math.Abs(by.y));
            bz = Vector3.Normalize(Vector3.Cross(refere, by));
            bx = Vector3.Normalize(Vector3.Cross(by, bz));
        }
    }
}
