// xvaProgramming.cpp : DLL �A�v���P�[�V�����p�ɃG�N�X�|�[�g�����֐����`���܂��B
//

#include "stdafx.h"
#include "xvaProgramming.h"


// ����́A�G�N�X�|�[�g���ꂽ�ϐ��̗�ł��B
XVAPROGRAMMING_API int nxvaProgramming=0;

// ����́A�G�N�X�|�[�g���ꂽ�֐��̗�ł��B
XVAPROGRAMMING_API int fnxvaProgramming(void)
{
    return 42;
}

// ����́A�G�N�X�|�[�g���ꂽ�N���X�̃R���X�g���N�^�[�ł��B
// �N���X��`�Ɋւ��Ă� xvaProgramming.h ���Q�Ƃ��Ă��������B
CxvaProgramming::CxvaProgramming()
{
    return;
}
