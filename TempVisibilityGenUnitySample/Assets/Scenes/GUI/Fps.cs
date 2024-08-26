using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine.Audio;
using UnityEngine.SceneManagement;
using UnityEngine.UI;
public class Fps : MonoBehaviour
{
    public Text fpsText;
    private float count;
    // public TextMeshProUGUI FpsText;

    private IEnumerator Start()
    {

        GUI.depth = 2;
        while (true)
        {
            count = 1f / Time.unscaledDeltaTime;
            yield return new WaitForSeconds(0.1f);
        }
    }

    private void OnGUI()
    {
        //GUI.Label(new Rect(5, 40, 100, 25), "FPS: " + Mathf.Round(count));
        fpsText.text = "Fps:" + Mathf.Round(count).ToString("#0");
    }
}