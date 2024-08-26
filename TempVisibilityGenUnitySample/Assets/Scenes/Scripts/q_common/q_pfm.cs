using System;
using System.IO;
using UnityEngine;
using static q_common.q_common;
using System.Text;

namespace q_common
{
    public class q_pfm
    {
        int w;
        int h;

        float[,] pixelsR;
        float[,] pixelsG;
        float[,] pixelsB;
        float[] colorStreamRGB;
        float[] colorStreamRGBA;

        public float[] ColorStreamRGB => colorStreamRGB;
        public float[] ColorStreamRGBA => colorStreamRGBA;

        public q_pfm() 
        {
            w = 0;
            h = 0;
        }
        public q_pfm(string filePath)
        {
            w = 0;
            h = 0;
            Load(filePath);
        }

        public void Load(string filePath)
        {
            using (var stream = new FileStream(filePath, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    // Read the file header
                    string magicNumber = new string(reader.ReadChars(2));
                    reader.ReadChar();
                    if (magicNumber != "PF" && magicNumber != "pf")
                    {
                        throw new Exception("Invalid PFM file format.");
                    }

                    string dimensionLine = ReadLine(reader);
                    string[] dimensions = dimensionLine.Trim().Split(' ');
                    w = int.Parse(dimensions[0]);
                    h = int.Parse(dimensions[1]);

                    // Check endianness
                    float scale = float.Parse(ReadLine(reader));
                    bool littleEndian = (scale < 0.0f);
                    // Calculate the absolute scale
                    scale = Math.Abs(scale);

                    // Read and parse pixel data
                    pixelsR = new float[w, h];
                    pixelsG = new float[w, h];
                    pixelsB = new float[w, h];

                    colorStreamRGB = new float[w * h * 3];
                    colorStreamRGBA = new float[w * h * 4];

                    for (int y = 0; y < h; y++) // Invert the image vertically
                    {
                        for (int x = 0; x < w; x++)
                        {
                            for (int z = 0; z < 3; z++)
                            {
                                float pixelValue = littleEndian
                                    ? reader.ReadSingle()
                                    : ReverseBytes(reader.ReadSingle());

                                if (z % 3 == 0)
                                    pixelsR[x, y] = pixelValue * scale;
                                else if (z % 3 == 1)
                                    pixelsG[x, y] = pixelValue * scale;
                                else if (z % 3 == 2)
                                    pixelsB[x, y] = pixelValue * scale;

                                colorStreamRGB[y * w + x + z] = pixelValue;
                            }
                        }
                    }

                    // Fill in colorStreamRGBA array
                    for (int y = 0; y < h; y++)
                    {
                        for (int x = 0; x < w; x++)
                        {
                            colorStreamRGBA[y * w + x] = pixelsR[x, y];
                            colorStreamRGBA[y * w + x + 1] = pixelsG[x, y];
                            colorStreamRGBA[y * w + x + 2] = pixelsB[x, y];
                            colorStreamRGBA[y * w + x + 3] = 1f;
                        }
                    }
                }
            }
        }

        public void Save(string filePath)
        {
            using (var stream = File.Open(filePath, FileMode.Create))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    string str = $"PF\n{(int)w} {(int)h}\n-1.000000\n";

                    byte[] bytes = Encoding.ASCII.GetBytes(str);
                    writer.Write(bytes);

                    for (int y = 0; y < h; y++) // Invert the image vertically
                    {
                        for (int x = 0; x < w; x++)
                        {
                            for (int z = 0; z < 3; z++)
                            {
                                if (z % 3 == 0)
                                    writer.Write(pixelsR[x, y]);
                                else if (z % 3 == 1)
                                    writer.Write(pixelsG[x, y]);
                                else if (z % 3 == 2)
                                    writer.Write(pixelsB[x, y]);
                            }
                        }
                    }
                }
                Debug.Log("The pfm is successfully saved to local: " + filePath);
            }
        }

    }
}
