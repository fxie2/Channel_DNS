template <typename T>
class vector
{
public:
	vector(int s = 0);
	~vector();
	//vector +-*/= vector
	vector<T>& operator+=(const vector<T> &v2);
	vector<T>& operator-=(const vector<T> &v2);
	vector<T>& operator*=(const vector<T> &v2);
	vector<T>& operator/=(const vector<T> &v2);
	friend vector<T> operator+(const vector<T> &v1, const vector<T> &v2);
	friend vector<T> operator-(const vector<T> &v1, const vector<T> &v2);
	friend vector<T> operator*(const vector<T> &v1, const vector<T> &v2);
	friend vector<T> operator/(const vector<T> &v1, const vector<T> &v2);
	//vector +-*/= scaler
	vector<T>& operator+=(double scaler);
	vector<T>& operator-=(double scaler);
	vector<T>& operator*=(double scaler);
	vector<T>& operator/=(double scaler);
	friend vector<T> operator+(const vector<T> &v, double scaler);
	friend vector<T> operator-(const vector<T> &v, double scaler);
	friend vector<T> operator*(const vector<T> &v, double scaler);
	friend vector<T> operator/(const vector<T> &v, double scaler);
	friend vector<T> operator+(double scaler, const vector<T> &v);
	friend vector<T> operator-(double scaler, const vector<T> &v);
	friend vector<T> operator*(double scaler, const vector<T> &v);
	friend vector<T> operator/(double scaler, const vector<T> &v);
	//Set value function
	vector<T>& operator=(const vector<T> &v2);
	vector<T>& operator=(double scaler);
	T& operator[](int rank);
private:
	T* _elem;
	int size;
};

template <typename T>
vector<T> operator+(vector<T> &v1, vector<T> &v2);

template <typename T>
vector<T> operator-(vector<T> &v1, vector<T> &v2);

template <typename T>
vector<T> operator*(vector<T> &v1, vector<T> &v2);

template <typename T>
vector<T> operator/(vector<T> &v1, vector<T> &v2);

template <typename T>
vector<T> operator+(const vector<T> &v, double scaler);

template <typename T>
vector<T> operator-(const vector<T> &v, double scaler);

template <typename T>
vector<T> operator*(const vector<T> &v, double scaler);

template <typename T>
vector<T> operator/(const vector<T> &v, double scaler);

template <typename T>
vector<T> operator+(double scaler, const vector<T> &v);

template <typename T>
vector<T> operator-(double scaler, const vector<T> &v);

template <typename T>
vector<T> operator*(double scaler, const vector<T> &v);

template <typename T>
vector<T> operator/(double scaler, const vector<T> &v);

