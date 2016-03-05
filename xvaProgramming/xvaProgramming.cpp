// xvaProgramming.cpp : DLL アプリケーション用にエクスポートされる関数を定義します。
//

#include "stdafx.h"
#include "xvaProgramming.h"


// これは、エクスポートされた変数の例です。
XVAPROGRAMMING_API int nxvaProgramming=0;

// これは、エクスポートされた関数の例です。
XVAPROGRAMMING_API int fnxvaProgramming(void)
{
    return 42;
}

// これは、エクスポートされたクラスのコンストラクターです。
// クラス定義に関しては xvaProgramming.h を参照してください。
CxvaProgramming::CxvaProgramming()
{
    return;
}
