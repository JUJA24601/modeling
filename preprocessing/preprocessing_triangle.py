import numpy as np
from stl import mesh

# stlファイルからcsvファイルを作成する

def main():
    triangles = mesh.Mesh.from_file("mesh.stl").vectors     # 三角形のデータ
    num_of_triangles = triangles.shape[0]        # 三角形の数
    print(triangles.shape)
    
    header = "apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z"
    csv = np.empty((num_of_triangles, 23))
    
    O = np.zeros(3)
    
    for i in range(num_of_triangles) :
        coord = triangles[i]
        side = np.array([distance(coord[0], coord[1]), distance(coord[1], coord[2]), distance(coord[2], coord[0])])
        a = side.max()
        p0 = coord[side.argmax()]
        b = distance(coord[(side.argmax()+2)%3], p0)
        # temp1 = coord[1] - coord[0]
        # temp2 = coord[2] - coord[0]
        # b = np.sqrt((distance(temp1,O)**2)*(distance(temp2,O)**2) - np.dot(temp1,temp2)**2)/a
        n1 = (coord[(side.argmax()+1)%3] - p0)/np.linalg.norm(coord[(side.argmax()+1)%3] - p0)
        n2 = (coord[(side.argmax()+2)%3] - p0)/np.linalg.norm(coord[(side.argmax()+2)%3] - p0)
        n3 = np.cross(n1,n2)/np.linalg.norm(np.cross(n1,n2))
        temp = np.concatenate([coord.ravel(), [a], [b], p0, n1, n2, n3])
        csv[i] = np.copy(temp)
        
    np.savetxt('data_triangle.csv', csv, delimiter=',', header=header, comments='')

def incenter(triangle) :
    """
    三角形の内心を求める
    
    Parameters
    ----------
    triangle : float (3,3)
        三角形の座標, triangle.shape -> (3,3) の必要あり
    
    Returns
    ---------
    incentre : float (1,3)
        三角形の内心
    """
    a = triangle[0]
    b = triangle[1]
    c = triangle[2]
    ab = distance(a,b)
    bc = distance(b,c)
    ca = distance(c,a)
    incentre = (bc*a+ca*b+ab*c)/(ab+bc+ca)
    return incentre

def distance(a, b) :
    """
    2点間 a, b の距離を求める
    
    Parameters
    ----------
    a : float (3,1)
        aの座標
    b : float (3,1)
        bの座標
    
    Returns
    ----------
    distance : float
        2点間の距離
    """
    distance = np.linalg.norm(b-a)
    
    return distance

if __name__ == '__main__':
    main()