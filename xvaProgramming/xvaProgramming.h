// 以下の ifdef ブロックは DLL からのエクスポートを容易にするマクロを作成するための 
// 一般的な方法です。この DLL 内のすべてのファイルは、コマンド ラインで定義された XVAPROGRAMMING_EXPORTS
// シンボルを使用してコンパイルされます。このシンボルは、この DLL を使用するプロジェクトでは定義できません。
// ソースファイルがこのファイルを含んでいる他のプロジェクトは、 
// XVAPROGRAMMING_API 関数を DLL からインポートされたと見なすのに対し、この DLL は、このマクロで定義された
// シンボルをエクスポートされたと見なします。
#ifdef XVAPROGRAMMING_EXPORTS
#define XVAPROGRAMMING_API __declspec(dllexport)
#else
#define XVAPROGRAMMING_API __declspec(dllimport)
#endif

// このクラスは xvaProgramming.dll からエクスポートされました。
class XVAPROGRAMMING_API CxvaProgramming {
public:
	CxvaProgramming(void);
	// TODO: メソッドをここに追加してください。
};

extern XVAPROGRAMMING_API int nxvaProgramming;

XVAPROGRAMMING_API int fnxvaProgramming(void);
