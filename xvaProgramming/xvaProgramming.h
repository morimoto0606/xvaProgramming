// �ȉ��� ifdef �u���b�N�� DLL ����̃G�N�X�|�[�g��e�Ղɂ���}�N�����쐬���邽�߂� 
// ��ʓI�ȕ��@�ł��B���� DLL ���̂��ׂẴt�@�C���́A�R�}���h ���C���Œ�`���ꂽ XVAPROGRAMMING_EXPORTS
// �V���{�����g�p���ăR���p�C������܂��B���̃V���{���́A���� DLL ���g�p����v���W�F�N�g�ł͒�`�ł��܂���B
// �\�[�X�t�@�C�������̃t�@�C�����܂�ł��鑼�̃v���W�F�N�g�́A 
// XVAPROGRAMMING_API �֐��� DLL ����C���|�[�g���ꂽ�ƌ��Ȃ��̂ɑ΂��A���� DLL �́A���̃}�N���Œ�`���ꂽ
// �V���{�����G�N�X�|�[�g���ꂽ�ƌ��Ȃ��܂��B
#ifdef XVAPROGRAMMING_EXPORTS
#define XVAPROGRAMMING_API __declspec(dllexport)
#else
#define XVAPROGRAMMING_API __declspec(dllimport)
#endif

// ���̃N���X�� xvaProgramming.dll ����G�N�X�|�[�g����܂����B
class XVAPROGRAMMING_API CxvaProgramming {
public:
	CxvaProgramming(void);
	// TODO: ���\�b�h�������ɒǉ����Ă��������B
};

extern XVAPROGRAMMING_API int nxvaProgramming;

XVAPROGRAMMING_API int fnxvaProgramming(void);
