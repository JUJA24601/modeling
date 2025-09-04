import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation
from matplotlib import animation, rc

def main():
    # # 円柱用の四角形を作成する
    # n = 64    # 何角柱
    # r = 0.5
    # r_poly = np.sqrt((2*np.pi*r**2)/(n*np.sin(2*np.pi/n)))
    # # print(f"r_poly : {r_poly}")
    
    # div = 100
    # rectangles = np.zeros(((2*div)*n,4,3))
    # print(rectangles.shape)
    # tik = np.zeros(2*div+1)
    # tik[0] = 0
    # tik[2*div] = 3
    # temp = 0
    # for i in range(div):
    #     # 端を細かくする場合
    #     r_0 = 1.5*(1 - (i/div)**2)
    #     r_1 = 1.5*(1 - ((i+1)/div)**2)

    #     # 均等に分割する場合
    #     # r_0 = 1.5*(1 - (i/div))
    #     # r_1 = 1.5*(1 - ((i+1)/div))

    #     temp += (r_0 - r_1)
    #     tik[2*div-1-i] = 3 - temp
    #     tik[i+1] = -(0 - temp)
        
    # count = 0

    # for i in range(2*div):
    #     x_1 = tik[i]
    #     x_2 = tik[i+1]
    #     # print(x_1, x_2)
    #     for j in range(n):
    #         y_1 = r_poly*np.sin(j*2*np.pi/n)
    #         y_2 = r_poly*np.sin((j+1)*2*np.pi/n)
    #         z_1 = r_poly*np.cos(j*2*np.pi/n)
    #         z_2 = r_poly*np.cos((j+1)*2*np.pi/n)
    #         # temp = np.array([[x_1,y_1,z_1],
    #         #                 [x_1,y_2,z_2],
    #         #                 [x_2,y_2,z_2],
    #         #                 [x_2,y_1,z_1]])
    #         temp = np.array([[y_1,z_1,x_1],
    #                          [y_2,z_2,x_1],
    #                          [y_2,z_2,x_2],
    #                          [y_1,z_1,x_2]])
    #         rectangles[count] = temp
    #         count += 1
    
    
    
    # # ここで四角形を作成する
    # div = 2
    # rectangles = np.zeros((6*((2*div)**2),4,3)) # こっちは立方体の場合
    # # rectangles = np.zeros((2*((2*div)**2),4,3)) # こっちは一枚の四角形の場合
    # print(rectangles.shape)
    # tik = np.zeros(2*div+1)
    
    # r = 2 # 四角形の半径
    
    # tik[0] = -r # 四角形のマイナス始点
    # tik[2*div] = r # 四角形のプラス終点
    # temp = 0
    # for i in range(div):
    #     # # 端を細かくする場合
    #     # r_0 = r*(1 - (i/div)**2)
    #     # r_1 = r*(1 - ((i+1)/div)**2)
        
    #     # 均等に分割する場合
    #     r_0 = r*(1 - (i/div))
    #     r_1 = r*(1 - ((i+1)/div))
        
    #     temp += (r_0 - r_1)
    #     tik[2*div-1-i] = r - temp
    #     tik[i+1] = -(r - temp)
        
    
    # count = 0
    # # face 1
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 64.5
    #         y_2 = tik[j+1] + 64.5
    #         temp = np.array([[x_1,y_1,-r+63],
    #                          [x_2,y_1,-r+63],
    #                          [x_2,y_2,-r+63],
    #                          [x_1,y_2,-r+63]])
    #         rectangles[count] = temp
    #         count += 1
            
    # # face 2
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 64.5
    #         y_2 = tik[j+1] + 64.5
    #         temp = np.array([[x_1,y_1,r+63],
    #                          [x_2,y_1,r+63],
    #                          [x_2,y_2,r+63],
    #                          [x_1,y_2,r+63]])
    #         rectangles[count] = temp
    #         count += 1
    
    # # face 3
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 63
    #         y_2 = tik[j+1] + 63
    #         temp = np.array([[x_1,-r+64.5,y_1],
    #                          [x_2,-r+64.5,y_1],
    #                          [x_2,-r+64.5,y_2],
    #                          [x_1,-r+64.5,y_2]])
    #         rectangles[count] = temp
    #         count += 1
            
    # # face 4
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 63
    #         y_2 = tik[j+1] + 63
    #         temp = np.array([[x_1,r+64.5,y_1],
    #                          [x_2,r+64.5,y_1],
    #                          [x_2,r+64.5,y_2],
    #                          [x_1,r+64.5,y_2]])
    #         rectangles[count] = temp
    #         count += 1
            
    # # face 5
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 63
    #         y_2 = tik[j+1] + 63
    #         temp = np.array([[-r+64.5,x_1,y_1],
    #                          [-r+64.5,x_2,y_1],
    #                          [-r+64.5,x_2,y_2],
    #                          [-r+64.5,x_1,y_2]])
    #         rectangles[count] = temp
    #         count += 1
            
    # # face 6
    # for i in range(2*div):
    #     x_1 = tik[i] + 64.5
    #     x_2 = tik[i+1] + 64.5
    #     # print(x_1, x_2)
    #     for j in range(2*div):
    #         y_1 = tik[j] + 63
    #         y_2 = tik[j+1] + 63
    #         temp = np.array([[r+64.5,x_1,y_1],
    #                          [r+64.5,x_2,y_1],
    #                          [r+64.5,x_2,y_2],
    #                          [r+64.5,x_1,y_2]])
    #         rectangles[count] = temp
    #         count += 1
            
    # num_of_rectangles = rectangles.shape[0]        # 四角形の数
    
    
    div = 8
    # rectangles = np.zeros((6*((2*div)**2),4,3))
    rectangles = np.zeros((((2*div)**2)*2+div*2*4,4,3))

    print(rectangles.shape)
    tik = np.zeros(2*div+1)
    r = 8
    tik[0] = -r
    tik[2*div] = r
    temp = 0

    for i in range(div):
        # 端を細かくする場合
        # r_0 = 0.5*(1 - (i/div)**2)
        # r_1 = 0.5*(1 - ((i+1)/div)**2)

        # 均等に分割する場合
        r_0 = r*(1 - (i/div))
        r_1 = r*(1 - ((i+1)/div))

        temp += (r_0 - r_1)
        tik[2*div-1-i] = r - temp
        tik[i+1] = -(r - temp)


    count = 0
    # face 1
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        for j in range(2*div):
            y_1 = tik[j] + 64.5
            y_2 = tik[j+1] + 64.5
            temp = np.array([[x_1,y_1,-1+63],
                            [x_2,y_1,-1+63],
                            [x_2,y_2,-1+63],
                            [x_1,y_2,-1+63]])
            rectangles[count] = temp
            count += 1

    # face 2
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        for j in range(2*div):
            y_1 = tik[j] + 64.5
            y_2 = tik[j+1] + 64.5
            temp = np.array([[x_1,y_1,0+63],
                            [x_2,y_1,0+63],
                            [x_2,y_2,0+63],
                            [x_1,y_2,0+63]])
            rectangles[count] = temp
            count += 1

    # face 3
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        y_1 = 62
        y_2 = 63
        temp = np.array([[x_1,-r+64.5,y_1],
                        [x_2,-r+64.5,y_1],
                        [x_2,-r+64.5,y_2],
                        [x_1,-r+64.5,y_2]])
        rectangles[count] = temp
        count += 1

    # face 4
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        y_1 = 62
        y_2 = 63
        temp = np.array([[x_1,r+64.5,y_1],
                        [x_2,r+64.5,y_1],
                        [x_2,r+64.5,y_2],
                        [x_1,r+64.5,y_2]])
        rectangles[count] = temp
        count += 1

    # face 5
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        y_1 = 62
        y_2 = 63
        temp = np.array([[-r+64.5,x_1,y_1],
                        [-r+64.5,x_2,y_1],
                        [-r+64.5,x_2,y_2],
                        [-r+64.5,x_1,y_2]])
        rectangles[count] = temp
        count += 1

    # face 6
    for i in range(2*div):
        x_1 = tik[i] + 64.5
        x_2 = tik[i+1] + 64.5
        # print(x_1, x_2)
        y_1 = 62
        y_2 = 63
        temp = np.array([[r+64.5,x_1,y_1],
                        [r+64.5,x_2,y_1],
                        [r+64.5,x_2,y_2],
                        [r+64.5,x_1,y_2]])
        rectangles[count] = temp
        count += 1

    num_of_rectangles = rectangles.shape[0]        # 4角形の数
    
    
    ################################################
    # 可視化
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Poly3DCollection を作成
    poly = Poly3DCollection(rectangles, facecolors='wheat', linewidths=0.1, edgecolors='saddlebrown', alpha=1.0)

    # Axes に追加
    ax.add_collection3d(poly)

    # 軸ラベルを設定
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim(-10+64.5, 10+64.5)
    ax.set_ylim(-10+64.5, 10+64.5)
    ax.set_zlim(-10+63, 10+63)
    ax.set_xticks([-10+64,-5+64, 0+64, 5+64, 10+64])
    ax.set_yticks([-10+64,-5+64, 0+64, 5+64, 10+64])
    ax.set_zticks([-10+63,-5+63, 0+63, 5+63, 10+63])
    ax.set_box_aspect([1, 1, 1])

    def update(frame):
        ax.view_init(elev=10, azim=frame)
        return poly,

    ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), blit=True)

    # ani.save('animation_rectangle.mp4', writer='ffmpeg', fps=15, dpi=350)
    ani.save('animation_rectangle.gif', fps=15, dpi=350)
    
    #################################################
    
    header = "apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, apex4_x, apex4_y, apex4_z, a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z"
    csv = np.empty((num_of_rectangles, 26))
    
    for i in range(num_of_rectangles) :
        coord = rectangles[i]
        a = distance(coord[0], coord[1])
        b = distance(coord[0], coord[3])
        p0 = coord[0]
        n1 = (coord[1] - coord[0])/a
        n2 = (coord[3] - coord[0])/b
        n3 = np.cross(n1,n2)
        temp = np.concatenate([coord.ravel(), [a], [b], p0, n1, n2, n3])
        csv[i] = np.copy(temp)
        
    np.savetxt('data_rectangle.csv', csv, delimiter=',', header=header, comments='')

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