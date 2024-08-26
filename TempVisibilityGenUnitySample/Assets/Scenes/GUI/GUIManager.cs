using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Linq;

#if UNITY_EDITOR
using UnityEditor;
#endif

public class GUIManager : MonoBehaviour
{
    [SerializeField] GenVisibilityVulkan vMono;
    [SerializeField] FreeCam fCam;
    [SerializeField] Slider shownVSlider;
    [SerializeField] Text shownVText;
    [SerializeField] Toggle enableDrawing;
    [SerializeField] Slider lineWSlider;
    [SerializeField] Slider ColorRSlider;
    [SerializeField] Slider ColorGSlider;
    [SerializeField] Slider ColorBSlider;

    [SerializeField] Dropdown AnimationDropdown;
    [SerializeField] Slider frameSlider;
    [SerializeField] Text frameText;
    [SerializeField] Button nextFBtn;
    [SerializeField] Button playBtn;
    [SerializeField] Button stopBtn;

    [SerializeField] Dropdown cameraDropdown;
    [SerializeField] Camera depthCam;
    [SerializeField] Camera freeCam;
    [SerializeField] GameObject resetBtn;

    bool hasUpdatedVis = false;

    void Start()
    {
        shownVSlider.value = vMono.shownVIdx;
        shownVText.text = vMono.shownVIdx.ToString();
        frameText.text = "";
        frameSlider.maxValue = 23;

        LineRenderer m_lineRenderer = vMono.lineRendererPrefab.GetComponent<LineRenderer>();
        ColorRSlider.value = m_lineRenderer.startColor.r;
        ColorGSlider.value = m_lineRenderer.startColor.g;
        ColorBSlider.value = m_lineRenderer.startColor.b;
        lineWSlider.value = 0.01f;// m_lineRenderer.startWidth;
    }

    public void OnDrawVisibility()
    {
        vMono.isDrawVisibility = enableDrawing.isOn;

        if (enableDrawing.isOn && !hasUpdatedVis)
        {
            vMono.OnShowFrame();
            hasUpdatedVis = true;
        }
        else
        { 
            vMono.DestroyVisibilityDraw();
            hasUpdatedVis = false;
        }
    }

    public void OnCameraChange()
    {
        if (cameraDropdown.value == 0)
        {
            depthCam.gameObject.SetActive(true);
            freeCam.gameObject.SetActive(false);
            resetBtn.gameObject.SetActive(false);
        }
        else
        {
            depthCam.gameObject.SetActive(false);
            freeCam.gameObject.SetActive(true);
            resetBtn.gameObject.SetActive(true);
        }
    }

    public void OnLineWidthChange()
    {
        LineRenderer m_lineRenderer = vMono.lineRendererPrefab.GetComponent<LineRenderer>();
        m_lineRenderer.startWidth = lineWSlider.value;
        m_lineRenderer.endWidth = lineWSlider.value;

        vMono.UpdateVisibilityWidth(lineWSlider.value);
    }

    public void OnColorRChange()
    {
        LineRenderer m_lineRenderer = vMono.lineRendererPrefab.GetComponent<LineRenderer>();

        Color newColor = new Color(
            ColorRSlider.value,
            m_lineRenderer.startColor.g, 
            m_lineRenderer.startColor.b);

        m_lineRenderer.startColor = newColor;
        m_lineRenderer.endColor = newColor;

        vMono.UpdateVisibilityColor(newColor);
    }

    public void OnColorGChange()
    {
        LineRenderer m_lineRenderer = vMono.lineRendererPrefab.GetComponent<LineRenderer>();

        Color newColor = new Color(
            m_lineRenderer.startColor.r,
            ColorGSlider.value,
            m_lineRenderer.startColor.b);

        m_lineRenderer.startColor = newColor;
        m_lineRenderer.endColor = newColor;

        vMono.UpdateVisibilityColor(newColor);
    }

    public void OnColorBChange()
    {
        LineRenderer m_lineRenderer = vMono.lineRendererPrefab.GetComponent<LineRenderer>();

        Color newColor = new Color(
            m_lineRenderer.startColor.r,
            m_lineRenderer.startColor.b,
            ColorBSlider.value);

        m_lineRenderer.startColor = newColor;
        m_lineRenderer.endColor = newColor;

        vMono.UpdateVisibilityColor(newColor);
    }


    public void ResetCameraView()
    {
        fCam.CameraReset();
    }

    public void OnShownVChange() 
    {
        vMono.shownVIdx = (int)shownVSlider.value;
        shownVText.text = ((int)shownVSlider.value).ToString();
    }

    private void Update()
    {
        if (vMono.isPlayingAnimation)
        {
            if (frameSlider.value == frameSlider.maxValue)
                frameSlider.value = 0;
            else
                frameSlider.value = vMono.currentFrame + 1;
            frameText.text = frameSlider.value.ToString();
        }
    }

    public void PlayAnimation()
    {
        vMono.isPlayingAnimation = true;

        if (frameSlider.value == frameSlider.maxValue)
            frameSlider.value = 0;
        else
            frameSlider.value = vMono.currentFrame + 1;
        frameText.text = frameSlider.value.ToString();
    }

    public void StopAnimation()
    {
        vMono.isPlayingAnimation = false;
    }

    public void OnAnimationIdxChange()
    {
        frameSlider.value = 0;
        frameText.text = "";

        vMono.currentFrame = -1;
        vMono.AnimIdx = AnimationDropdown.value;
        vMono.ResetModel();

        frameSlider.maxValue = vMono.frameDatasWS[vMono.AnimIdx].Length;
    }

    public void NextFrameAnimation()
    {
        vMono.playNextFrame = true;
        if (frameSlider.value == frameSlider.maxValue - 1)
            frameSlider.value = 0;
        else
            frameSlider.value = vMono.currentFrame + 1;
        frameText.text = frameSlider.value.ToString();
    }

    public void ExitGame()
    {
#if UNITY_EDITOR
        UnityEditor.EditorApplication.isPlaying = false;
#else
        Application.Quit();
#endif
    }

}
