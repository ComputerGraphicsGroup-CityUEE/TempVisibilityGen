using UnityEngine;
using UnityEngine.EventSystems;

public enum CamViewMode
{
    FREE,
    TOP,
    LIMIT
}
public class FreeCam : MonoBehaviour
{
    private Transform camTrans;
    private bool _isFirstClick = true;


    CamViewMode viewMode = CamViewMode.FREE;

    [SerializeField]
    private Vector3 _resetTrans; // Camera Reset Position
    [SerializeField]
    private Vector3 _resetAngles;// Camera Reset Rotation

    [Header("Middle mouse button movement speed")]
    public float m_mSpeed = 0.5f;
    [Header("Rotate Speed")]
    public float m_rSpeed = 5f;
    [Header("Scale Speed")]
    public float m_sSpeed = 5f;
    [Header("Max Scale Distance")]
    public float m_maxDistance = 10f;
    [Header("Easing value for middle button movement")]
    public float moveSmoothing = 0.2f;

    private float m_deltX = 0f;
    private float m_deltY = 0f;

    void Start()
    {
        camTrans = transform;
    }
    void Update()
    {
        if (EventSystem.current.IsPointerOverGameObject()) return;

        if (viewMode != CamViewMode.LIMIT)
        {
            CameraMouseControl();
        }
    }
    private void CameraMouseControl()
    {
        // camera scale
        if (Input.GetAxis("Mouse ScrollWheel") != 0)
        {
            float m_distance = Input.GetAxis("Mouse ScrollWheel") * m_sSpeed;
            Vector3 newPos = camTrans.localPosition + camTrans.forward * m_distance;
            if (newPos.magnitude >= m_maxDistance) return;
            camTrans.localPosition = newPos;
        }
        // Camera rotation;
        else if (Input.GetMouseButton(1))
        {
            if (!_isFirstClick)
            {
                m_deltX += Input.GetAxis("Mouse X") * m_rSpeed;
                m_deltY -= Input.GetAxis("Mouse Y") * m_rSpeed;
            }
            else 
            {
                _isFirstClick = false;
                m_deltX = _resetAngles.y;
                m_deltY = _resetAngles.x;
            }

            m_deltX = ClampAngle(m_deltX, -360, 360);
            m_deltY = ClampAngle(m_deltY, -70, 70);

            camTrans.localRotation = Quaternion.Euler(m_deltY, m_deltX, 0);
        }
        // camera move
        else if (Input.GetMouseButton(2))
        {
            float dx = Input.GetAxis("Mouse X");
            float dy = Input.GetAxis("Mouse Y");
            Vector3 yz = camTrans.forward + camTrans.up;
            yz.y = 0;
            Vector3 TargetLookAt = camTrans.position;
            TargetLookAt -= (yz * dy + transform.right * dx) * m_mSpeed;
            camTrans.position = Vector3.Lerp(camTrans.position, TargetLookAt, moveSmoothing);
        }
    }

    public void CameraReset()
    {
        camTrans.localPosition = _resetTrans;
        camTrans.localRotation = Quaternion.Euler(_resetAngles);

    }

    float ClampAngle(float angle, float minAngle, float maxAgnle)
    {
        if (angle <= -360)
            angle += 360;
        if (angle >= 360)
            angle -= 360;

        return Mathf.Clamp(angle, minAngle, maxAgnle);
    }

}